import time
import numpy as np
from collections import defaultdict
import itertools
import copy
import random
import pandas as pd
import memory_profiler
import math

alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
            'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
df = pd.read_csv("Diabetes-Data/data-01.csv",
                 sep='\t', header=None, usecols=[2])

# reading multiple files for diabetes. Ignore this if you want to use your own dataset
for i in range(1, 20):
    if i < 10:
        df2 = pd.read_csv("Diabetes-Data/data-0"+str(i)+".csv",
                          sep='\t', header=None, usecols=[2])
    else:
        df2 = pd.read_csv("Diabetes-Data/data-"+str(i)+".csv",
                          sep='\t', header=None, usecols=[2])
    df = pd.concat([df, df2], ignore_index=True)


itemssss = {}
jj = 0
qq = 0
S = "" # holds the entire string of the data
for row in df.iterrows():
    if qq == 600: # set the number of rows to use for mining patterns
        break
    qq += 1
    if str(int(row[1].values)) in itemssss:
        S += itemssss[str(int(row[1].values))]
    else:
        itemssss[str(int(row[1].values))] = alphabet[jj]
        S += itemssss[str(int(row[1].values))]
        jj += 1


def slice_events(events, keys):
    '''Remove specific keys from a dictionary'''
    for item in keys:
        del events[item]
    return events


def find_periodicity(event, occ_vec, k, threshhold):
    '''Checks if an event is periodic for a specific time'''
    op = []
    for i in range(1, k):
        st = occ_vec[i-1]
        for j in range(i, k):
            d = occ_vec[j]
            period = d - st
            if period == 0:
                continue
            v = []
            v.append(st)
            v.append(d)
            for m in range(j+1, k):
                if(occ_vec[m]-st) % period == 0:
                    v.append(occ_vec[m])
            pp = math.floor(((len(S)-st)+1)/period)
            if pp == 0:
                continue
            else:
                conf = len(v)/pp
                if conf > threshhold:
                    op.extend(v)
    return list(set(op))


def calc_pattern2(single_events, double_events, prev_keys, max_pattern_len):
    '''Calculates patterns for the newly discoverd pattern'''
    new_pattern = defaultdict(list)
    new_events = copy.deepcopy(single_events)
    for k, v in double_events.items():
        periodic = find_periodicity(k, v, len(v), 0.6)
        for i in range(len(periodic)):
            for key in new_events.keys():
                for items in new_events[key]:
                    diff = abs(periodic[i] - items)
                    if diff >= max_pattern_len or diff < len(k):
                        continue
                    if items > periodic[i]: # first event occurs before the second event. E.g for event A and B, this will create patterns starting with A
                        if diff > len(k):
                            star_count = diff-len(k)
                            if star_count > 3:
                                continue
                            stars = ""
                            for j in range(star_count):
                                stars += "*"
                            if periodic[i] in new_pattern[k+stars+key]:
                                continue
                            else:
                                if key != k:
                                    new_pattern[k+stars +
                                                key].append(periodic[i])
                        else:
                            if key != k:
                                new_pattern[k +
                                            key].append(periodic[i])
                    elif periodic[i] > items: # second event occurs before the first event. E.g for event A and B, this will create patterns starting with B
                        if diff > len(key):
                            star_count = diff-len(key)
                            if star_count > 3:
                                continue
                            stars = ""
                            for j in range(star_count):
                                stars += "*"
                            if periodic[i] in new_pattern[key+stars+k]:
                                continue
                            else:
                                if key != k:
                                    new_pattern[key+stars +
                                                k].append(items)
                        else:
                            if key != k:
                                new_pattern[key+k].append(items)
    return new_pattern


def calc_pattern(event1, occ_vec1, events, max_pattern_len, prev_keys, p):
    '''Calculates patterns for periodic events'''
    new_pattern = defaultdict(list)
    new_events = copy.deepcopy(events)
    new_events = slice_events(new_events, prev_keys) # remove the keys that we have previously mined patterns for
    looped = []
    periodic = find_periodicity(event1, occ_vec1, len(occ_vec1), 0.6) # check periodicity of the event
    for i in range(len(periodic)): # loop for every period that the event is perodic for
        for key in new_events.keys():
            if key != event1 or key == event1:
                for items in new_events[key]:
                    diff = abs(periodic[i] - items) # get the difference of occurence vectors
                    if diff >= max_pattern_len: # check against max pattern length
                        continue
                    if items > periodic[i]: # first event occurs before the second event. E.g for event A and B, this will create patterns starting with A
                        if diff > len(key):
                            star_count = diff-len(key)
                            if star_count > 3:
                                continue
                            stars = ""
                            for j in range(star_count):
                                stars += "*"
                            if periodic[i] not in new_pattern[event1+stars+key]:
                                new_pattern[event1+stars +
                                            key].append(periodic[i])
                        elif diff == 0:
                            new_pattern[event1 +
                                        key].append(periodic[i])
                        else:
                            if periodic[i] not in new_pattern[event1+key]:
                                new_pattern[event1 +
                                            key].append(periodic[i])
                    elif periodic[i] > items: # second event occurs before the first event. E.g for event A and B, this will create patterns starting with B
                        if diff > len(key):
                            star_count = diff-len(key)
                            if star_count > 3:
                                continue
                            stars = ""
                            for j in range(star_count):
                                stars += "*"
                            if items not in new_pattern[key+stars+event1]:
                                new_pattern[key+stars+event1].append(items)
                        elif diff == 0:
                            new_pattern[key+event1].append(items)
                        else:
                            new_pattern[key+event1].append(items)
    for z in range(p+1, len(new_events)): # calculate patterns for the newly found pattern
        new_pattern.update(calc_pattern2(
            new_events, new_pattern, looped, max_pattern_len))
    return new_pattern

#@profile
def mainfunc():
    prev_keys = []
    start = time.perf_counter() # initializing the calculation of time required to mine patterns 
    items = defaultdict(list)
    for pos, item in enumerate(S):
        items[item].append(pos)
        #print(items)

    patterns = {}
    p = 0
    for key, value in items.items():
        new_patterns = calc_pattern(key, value, items, 6, prev_keys, p)
        prev_keys.append(key)
        p += 1

    finish = time.perf_counter()
    print(finish-start)


if __name__ == '__main__':
    mainfunc()
