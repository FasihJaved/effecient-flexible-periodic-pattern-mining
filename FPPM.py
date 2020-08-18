import pandas as pd
import time
import numpy as np
from collections import defaultdict
import itertools
import copy
import random
import memory_profiler
import math

# mapping number ranges to alphabets
alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
            'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
            
df = pd.read_csv("hour.csv", sep=',', usecols=[16])
df = df.reindex(columns = ['cnt', 'cat'])  

for index, row in df.iterrows():
    if row['cnt'] <= 10:
        df['cat'][index] = 'a'
        continue
    if row['cnt'] <= 20:
        df['cat'][index] = 'b'
        continue
    if row['cnt'] <= 30:
        df['cat'][index] = 'c'
        continue
    if row['cnt'] <= 40:
        df['cat'][index] = 'd'
        continue
    if row['cnt'] <= 50:
        df['cat'][index] = 'e'
        continue
    if row['cnt'] <= 60:
        df['cat'][index] = 'f'
        continue
    if row['cnt'] <= 70:
        df['cat'][index] = 'g'
        continue
    if row['cnt'] <= 80:
        df['cat'][index] = 'h'
        continue
    if row['cnt'] <= 90:
        df['cat'][index] = 'i'
        continue
    if row['cnt'] <= 100:
        df['cat'][index] = 'j'
        continue
    if row['cnt'] <= 110:
        df['cat'][index] = 'k'
        continue
    if row['cnt'] <= 120:
        df['cat'][index] = 'l'
        continue
    if row['cnt'] <= 130:
        df['cat'][index] = 'm'
        continue
    if row['cnt'] <= 140:
        df['cat'][index] = 'n'
        continue
    if row['cnt'] <= 150:
        df['cat'][index] = 'o'
        continue
    if row['cnt'] <= 160:
        df['cat'][index] = 'p'
        continue
    if row['cnt'] <= 170:
        df['cat'][index] = 'q'
        continue
    if row['cnt'] <= 180:
        df['cat'][index] = 'r'
        continue
    if row['cnt'] <= 190:
        df['cat'][index] = 's'
        continue
    if row['cnt'] <= 200:
        df['cat'][index] = 't'
        continue
    if row['cnt'] <= 210:
        df['cat'][index] = 'u'
        continue
    if row['cnt'] <= 220:
        df['cat'][index] = 'v'
        continue
    if row['cnt'] <= 230:
        df['cat'][index] = 'w'
        continue
    if row['cnt'] <= 240:
        df['cat'][index] = 'x'
        continue
    if row['cnt'] <= 250:
        df['cat'][index] = 'y'
        continue
    if row['cnt'] <= 260:
        df['cat'][index] = 'z'
        continue

itemssss = {}
jj = 0
qq = 0
S = ""
for index, row in df.iterrows():
    if qq == 2000: # set the number of rows to mine e.g 2000 = first 2000 rows of a csv file
        break
    qq += 1
    S += str(row['cat'])


class SuffixTrie(object):
    def __init__(self, t):
        t += '$'  # add terminator
        self.root = {}
        self.occ_vec = {}
        self.len_vec = {}
        for i in range(len(t)):
            cur = self.root
            occ_vec = self.occ_vec
            already_passed = False
            inside_loop = i-1
            for c in t[i:]:
                inside_loop += 1
                if c not in cur:
                    cur[c] = {}
                if c not in occ_vec and already_passed == False:
                    occ_vec[c] = []
                    occ_vec[c].append(i)
                    self.len_vec[c] = len(t)-i
                    already_passed = True
                elif c in occ_vec:
                    added = False
                    for item in occ_vec[c]:
                        if item >= inside_loop:
                            added = True
                    if not added:
                        occ_vec[c].append(inside_loop)
                cur = cur[c]
                already_passed = True

    def followPath(self, s):
        '''Traverses the tree'''
        cur = self.root
        for c in s:
            if c not in cur:
                return None
            cur = cur[c]
            print(cur)
        return cur

    def hasSubstring(self, s):
        return self.followPath(s) is not None

    def hasSuffix(self, s):
        node = self.followPath(s)
        return node is not None and '$' in node


def find_periodicity(event, occ_vec, k, threshhold):
    '''Checks if an event is periodic for a specific time'''
    op = {}
    periodic = False
    for i in range(1, k):
        st = occ_vec[i-1]
        for j in range(i, k):
            d = occ_vec[j]
            period = d - st
            v = []
            v.append(st)
            v.append(d)
            for m in range(j+1, k):
                if(occ_vec[m]-st) % period == 0:
                    v.append(occ_vec[m])
            pp = math.floor((len(S)-st)/period)
            conf = len(v)/pp
            if conf > threshhold:
                periodic = True
    return periodic



def joianable(occ_vec1, occ_vec2, dist):
    '''Checks if occurence vectors have a higher distance than the max length of unimportant events allowed'''
    new_occ_vec = []
    for item in occ_vec1:
        for item2 in occ_vec2:
            if item2 == (item+dist):
                new_occ_vec.append(item)
    if len(new_occ_vec) > 0:
        return new_occ_vec
    else:
        return 0


def recursive_items(dictionary):
    '''Yeilds the dictionary items'''
    for key, value in dictionary.items():
        if type(value) is dict:
            yield (key, value)
            # yield from recursive_items(value)
        else:
            yield (key, value)


def trailingstarcount(word):
    '''Adds a star in pattern'''
    stars = 0
    for char in word:
        if char == '*':
            stars += 1

    return stars


def lad_fact(occ_vec, conf):
    '''Calculates the ladder factor'''
    ladder_factor = math.ceil(len(occ_vec)*conf)
    return abs(occ_vec[ladder_factor]-len(S))


def calculations(key, value, s):
    st = find_periodicity(key, s.occ_vec[key], len(s.occ_vec[key]), 0.60) # implementation of periodicity caulcualtion algo from base paper
    if st is False:
        return
    pt = {}
    current = []
    current.append(key)
    current_node = {key: value}
    l = 25  # ladder factor, determines the number of nodes to traverse per branch
    level = 0
    while l > 0: # main loop per level
        level += 1
        next_node = {}
        current_level = []
        for key2, value2 in current_node.items():
            next_node.update(value2.items())
            for keys, items in value2.items(): # checking children's periodicity
                st2 = find_periodicity(
                    keys, s.occ_vec[keys], len(s.occ_vec[keys]), 0.60)
                if st2 is True:
                    current_level.append(keys)
        Next = []
        for item in current_level: # setting for next level
            for item2 in current[level-1:]:
                if item != "$": # checking for root node
                    # checks the joinability of items based on their occurence vectors
                    if len(item) > 1 and len(item2) == 1:
                        joinable = joianable(
                            s.occ_vec[item2], pt[item], len(item))
                    elif len(item) == 1 and len(item2) > 1:
                        joinable = joianable(
                            pt[item2], s.occ_vec[item], len(item2))
                    elif len(item) > 1 and len(item2) > 1:
                        joinable = joianable(
                            pt[item2], pt[item], len(item2))
                    else:
                        joinable = joianable(
                            s.occ_vec[item2], s.occ_vec[item], len(item2))

                    # checking the periodicity of the new pattern
                    if type(joinable) == list:
                        pt[item2+item] = joinable
                        st2 = find_periodicity(
                            item2+item, joinable, len(joinable), 0.60)
                        if st2 is False:
                            del pt[item2+item]
                    
                    # checking which pattern to have stars
                    if l > 1:
                        if len(item2) > 1:
                            pt2 = {item2+"*": pt[item2]}
                        else:
                            pt2 = {item2+"*": s.occ_vec[item2]}
                        if pt.get(item2+"*") is None and trailingstarcount(item2+"*") < 3:
                            pt[item2+"*"] = pt2[item2+"*"]
        l -= 1
        current_node = next_node
        for item in pt:
            current.append(item)
    return pt

#@profile
def mainfunc():
    start = time.perf_counter() # initializing the calculation of time required to mine patterns 
    s = SuffixTrie(S)
    final_ptrns = []
    ptrn_index = []
    for key, value in s.root.items():
        pp = calculations(key, value, s) # calculations of occurence vector etc
        if pp is not None:
            for k, v in pp.items():
                if len(v)>0:
                    final_ptrns.append(k)
                    ptrn_index.append(v)
    
    #save file if needed
    #df = pd.DataFrame(list(zip(final_ptrns, ptrn_index)), columns=['Patterns', 'Indexes'])
    #df.to_csv("finalpatrns_base_2000_2.csv")

    finish = time.perf_counter()
    print(finish-start)


if __name__ == '__main__':
    mainfunc()

