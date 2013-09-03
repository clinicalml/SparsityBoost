'''
Created on Apr 29, 2013

@author: eliotpbrenner
'''
import bisect


def find_ge(a, x):
    'Find leftmost item, and its index, greater than or equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return a[i],i
    raise ValueError

def find_gt(a, x):
    'Find leftmost value, and its index, greater than x'
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return a[i],i
    raise ValueError

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x)
    if i:
        return a[i-1], i
    raise ValueError

def find_le_withBounds(a,x,lo,hi):
    'Find rightmost value less than or equal to x'
    'Considers subset of the list a, a[lo,hi]'
    i = bisect.bisect_right(a, x, lo, hi)
    if i:
        return a[i-1],i
    raise ValueError


def firstIndexHash(entryList, increasingList):
    #:Input:
	#  entryList: a list of entries from which the increasing list is drawn
	#  increasingList: an increasing list of entries, taken from entryList
	#:Output:
	#  Hash mapping each element of entryList to the leftmost index of
    #  increasingList where it appears
	return { entry : find_ge(increasingList, entry)[1] for entry in entryList }
