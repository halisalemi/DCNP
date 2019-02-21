#ifndef BINARYHEAP_H
#define BINARYHEAP_H
#include<iostream>
#include<climits>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// Structure to represent a binary heap node
struct BinaryHeapNode
{
	long  v;
	long dist;
};

// Structure to represent a binary heap
struct BinaryHeap
{
	long size;      // Number of heap nodes present currently
	long capacity;  // Capacity of min heap
	long *pos;     // This is needed for decreaseKey()
	struct BinaryHeapNode **array;
};


// A function to create a new Min Heap Node
struct BinaryHeapNode* newMinHeapNode(long v, long dist);

// A function to create a Min Heap
struct BinaryHeap* createMinHeap(long capacity);

// A function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct BinaryHeapNode** a, struct BinaryHeapNode** b);


// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct BinaryHeap* minHeap, long idx);


// A function to check if the given minHeap is ampty or not
long isEmpty(struct BinaryHeap* minHeap);


// Standard function to extract minimum node from heap
struct BinaryHeapNode* extractMin(struct BinaryHeap* minHeap);

// A function to delete minimum node from heap
void deleteMin(struct BinaryHeap* minHeap);


// Function to decreasy dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct BinaryHeap* minHeap, long v, long dist);


// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct BinaryHeap *minHeap, long v);



#endif

