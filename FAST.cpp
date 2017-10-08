#include "FAST.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

#include "prng.h"
/**********************************************************************************************
Implementation of the FAST algorithm to Find Frequent Items on weighted streams
Based on the paper Fast Flow Volume Estimation (ACM ICDCN 2018) 
of Ben Basat, Einziger and Friedman, 2016
Implementation by Ben Basat and Einziger, 2016-2017

Original Code: 2016-06
This version: 2017-10

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA.
***********************************************************************************************/

FAST::FAST(float fPhi, int M, float gamma)
{
	/*
		FAST's constructure, takes three params:
			1. fPhi - the error guarantee. The algorithm always has an error <= N*M*fPhi
			2. M - the maximal possible packet size
			3. Gamma - the additional space parameter. Higher Gamma => faster but uses more space
	*/
	int i;
	m_nrCounters = ceil((int) 1.0 / fPhi);
	m_M = M;
	m_gamma = gamma;
	m_stepSize = 1 + M * gamma / 2;
	m_stepSizeMinusOne = m_stepSize - 1;
	assert(m_stepSize);
	assert(gamma > 0);
	m_nrCounters *= (1 + gamma);

	m_a = (long long)698124007;
	m_b = (long long)5125833;
	assert(m_nrCounters > 0);
#if UNDER_ESTIMATOR
	m_maxRootVictim = 0;
#else
	m_maxRootVictim = m_stepSizeMinusOne;
#endif

	m_tblsz = FAST_HASHMULT*m_nrCounters;
	m_hashtable = (FASTITEM**)calloc(m_tblsz, sizeof(FASTITEM *));
	m_groups = (FASTGROUP *)calloc(m_nrCounters, sizeof(FASTGROUP));
	m_items = (FASTITEM *)calloc(m_nrCounters, sizeof(FASTITEM));
	m_freegroups = (FASTGROUP **)calloc(m_nrCounters, sizeof(FASTGROUP*));

	for (i = 0; i<m_tblsz; i++)
		m_hashtable[i] = NULL;

	m_root = m_groups;
	m_groups->count = 0;
	m_groups->nextg = NULL;
	m_groups->previousg = NULL;

	m_groups->items = m_items;
	for (i = 0; i<m_nrCounters; i++)
		m_freegroups[i] = &m_groups[i];
	m_gpt = 1; // initialize list of free groups

	for (i = 0; i<m_nrCounters; i++)
	{
		m_items[i].item = 0;
		m_items[i].delta = -1;
		m_items[i].hash = 0;
		m_items[i].nexti = NULL;
		m_items[i].previousi = NULL;  // initialize values

		m_items[i].parentg = m_groups;
		m_items[i].nexting = &(m_items[i + 1]);
		m_items[i].previousing = &(m_items[i - 1]);
		// create doubly linked list
	}
	m_items[0].previousing = &(m_items[m_nrCounters - 1]);
	m_items[m_nrCounters - 1].nexting = &(m_items[0]);
	// fix start and end of linked list

	m_nrFree = m_nrCounters;
}

FAST::~FAST()
{
	free(m_items);
	free(m_groups);
	free(m_freegroups);
	free(m_hashtable);
}



void FAST::update(int newitem, int weight) {
	// The public `update` function. Adds weight (weight) to flow (newitem)
	int h = getHash(newitem);
	FASTITEM *til = getCounter(newitem, h); //returns a pointer to newitem's counter if such exists or NULL otherwise
	FASTITEM &il = (til ? *til : takoverMinimal(newitem, h)); // gets a reference to the counter. Takes over a counter if needed.
#if UNDER_ESTIMATOR //Can be used to reduce the error
	if (!(til || m_nrFree))
			il.delta = m_root->count;
#endif
	int newVal = weight + (til ? il.remainder : m_maxRootVictim); // The newVal is set to the previous remainder + the new weight
	unsigned int nrIncs = newVal / m_stepSize; // How many groups in the SOS we need advance the current counter

	il.remainder = newVal % m_stepSize; //TODO: consider changing % to shifts
	if (nrIncs && (!m_nrFree || til)) { 
		//This is the path usually taken. We end here unless the update can be resolved just by updating the remainder 
		//or by allocating a free counter.
		advanceCounter(il, nrIncs);
		return;
	}
	if (m_nrFree && !til)
	{ 	// We deliberately use duplicated code as this `else' is rarely accessed. It is only accessed if newitem has no counter
		// and we have not yet allocated all counters.
		newVal = weight;
		nrIncs = newVal / m_stepSize;
		--m_nrFree;
		il.parentg->items = il.parentg->items->nexting;
		il.remainder = newVal % m_stepSize; //TODO: consider changing % to shifts
		if (nrIncs)
			advanceCounter(il, nrIncs);
	}
}

inline int const FAST::getHash(unsigned int item) {
	return hash31(m_a, m_b, item) % m_tblsz;
}

inline FASTITEM * FAST::getCounter(unsigned int item, int h) {
	// Returns a pointer to item's counter or NULL if such does not exists
	FASTITEM * il = m_hashtable[h];
	while (il && (il->item != item))
		il = il->nexti;
	return il;
}

inline FASTITEM & FAST::takoverMinimal(unsigned int item, int h) {
	// Takes a counter from the minimal-value group in SOS and reallocates it to `item`.
	FASTITEM & il = *(m_root->items);
	removeFromHash(il); // Remove the association with the item previously allocated to that counter.
	insertIntoHashtable(il, h, item);
#if UNDER_ESTIMATOR
		if (m_maxRootVictim < il.remainder)
		{
			m_maxRootVictim = il.remainder;
		}
#endif
	return il;
}

void FAST::advanceCounter(FASTITEM &il, unsigned int nrIncs) {
	// The function takes two parameters: a reference to a counter (il) and the number of increments we need to make (nrIncs)
	FASTGROUP &oldgroup = *(il.parentg); //The SOS group that the counter belonged to
	unsigned int goalVal = oldgroup.count + nrIncs; //The value that the counter should have at the end of the process.
	FASTGROUP &group = getLastGroup(oldgroup, goalVal); //The last group with a value <= goalVal

	if ((il.nexting) == &il) {      // if the counter is the only one in its group

		if (&group == &oldgroup) {    // if we can simply increase the oldgroup's count
			oldgroup.count = goalVal;
			return;
		}
		if (group.count == goalVal) { // if there exists a group with count = goalVal
			putInNewGroup(il, group);
			return;
		}
		moveGroup(il, group, goalVal);
		return;
	}
	disconnectCounter(il); // Remove the association from the old group
	assert(group.count <= goalVal);

	if (group.count == goalVal) {
		putInNewGroup(il, group); // if the goal group exists	
		return;
	}
	newGroup(il, group, goalVal); // Create a new group as no other counter share the same value
}

inline void FAST::disconnectCounter(FASTITEM & il) { 
	//Disconnects only if there exists other counters in the same group
	if (il.parentg->items == &il) {
		il.parentg->items = il.nexting;
	}
	assert(il.previousing);
	assert(il.nexting);
	il.previousing->nexting = il.nexting;
	il.nexting->previousing = il.previousing;
	il.nexting = &il;
	il.previousing = &il;
}

inline FASTGROUP & FAST::getLastGroup(FASTGROUP & oldgroup, unsigned int goalVal) {
	// Returns a reference to the last SOS group with a value <= goalVal
	FASTGROUP *group = &oldgroup;
	while (group->nextg && (group->nextg->count <= goalVal)) {
		group = group->nextg;
	}
	assert((group->items->nexting == group->items) == (group->items->previousing == group->items));
	return *group;
}

inline void FAST::removeFromHash(FASTITEM & il) {
	//Our hash table is chain-based and this function removes the item from the chain.
	
	// if il was first - the chain will point on the next item.
	if (m_hashtable[il.hash] == &il)
		m_hashtable[il.hash] = il.nexti;
	// if there is another node with the same hash. 
	if (il.nexti)
	{
		// next item - previous - points to my previous.
		il.nexti->previousi = il.previousi;

	}
	if (il.previousi) {
		//prev item next - points to my next.
		il.previousi->nexti = il.nexti;
	}
}

/*
inline void FAST::addToGroup(FASTGROUP & group, FASTITEM & il) {
	il.previousing = group.items->previousing;
	il.nexting = group.items;
	il.parentg = &group;
	group.items->previousing->nexti = &il;
	group.items->previousing = &il;
}*/



inline void const FAST::showGroups() {
	// A debug function that outputs the current SOS structure.
	FASTGROUP *g;
	FASTITEM *i, *first;
	int n, wt;

	g = m_groups;
	wt = 0;
	n = 0;
	while (g != NULL)
	{
		printf("Group %lld :", g->count);
		first = g->items;
		i = first;
		if (i != NULL)
			do
			{
				printf("%d -> ", i->item);
				i = i->nexting;
				wt += g->count;
				n++;
			} while (i != first);
		else printf(" empty");
		printf(")");
		g = g->nextg;
		if ((g != NULL) && (g->previousg->nextg != g))
			printf("Badly linked");
		printf("\n");
	}
	printf("In total, %d items, with a total count of %d\n", n, wt);
}

void FAST::insertIntoHashtable(FASTITEM &newi, int i, unsigned int newitem) {
	// Add counter (newi) with hash value i and ID newitem into the hash
	newi.nexti = m_hashtable[i];
	newi.item = newitem; // overwrite the old item
	newi.hash = i;
	newi.previousi = NULL;
	// insert item into the hashtable
	if (m_hashtable[i])
		m_hashtable[i]->previousi = &newi;
	m_hashtable[i] = &newi;
}

unsigned int const FAST::query(unsigned int x) {
	// The main query function. Given an identifier `x`, this returns its estimated volume.
	int h;
	FASTITEM *il;
	h = hash31(m_a, m_b, x) % m_tblsz;
	il = m_hashtable[h];
	// if x has a counter.
	while (il && il->item != x)
		il = il->nexti;
	if (il){
		if (il->delta == -1)
			return il->parentg->count * m_stepSize + il->remainder;
		return il->parentg->count * m_stepSize + il->remainder - (il->delta + 1)*m_stepSize - 1;
	}
	if (UNDER_ESTIMATOR)
		return 0;
	int minCount = m_root->count;
	// if all counters are used. 
	if (minCount)
		return minCount * m_stepSize - 1;
	// if there are unused counters and x does not have a counter, we are certain that x never arrived. 
	return 0;
}


inline void FAST::connectToGroup(FASTITEM &newi, FASTGROUP &tmpg) {
	// Associate the counter(newi) with the group (tmpg)
	newi.nexting = tmpg.items;
	newi.previousing = tmpg.items->previousing;
	assert((newi.previousing != &newi) && (newi.nexting != &newi));
	newi.previousing->nexting = &newi;
	newi.nexting->previousing = &newi;
}

inline void FAST::putInNewGroup(FASTITEM &newi, FASTGROUP & tmpg) {
	// Connects the counter (newi) to the group (tmpg) and recycles its previous group if it is now empty
	FASTGROUP * oldgroup = newi.parentg;
	assert(oldgroup != &tmpg);
	// put item in the tmpg group
	newi.parentg = &tmpg;
	assert((newi.previousing == &newi) && (newi.nexting == &newi));
	if (oldgroup->items == &newi) { // group will be empty
		recycleGroup(*oldgroup);
	}
	connectToGroup(newi, tmpg);
}

inline void FAST::moveGroup(FASTITEM &il, FASTGROUP &group, unsigned int goalVal) {
	// Used in the case where we need to move the group (group) for setting the value goalVal and keeping the SOS sorted
	assert(il.parentg->nextg && il.parentg->nextg->count < goalVal);
	FASTGROUP &tmpg = *(il.parentg);

	if (&tmpg == m_root) // May fail if we allocate only a single counter
	{
		m_root = tmpg.nextg;
		m_root->previousg = NULL;
#if UNDER_ESTIMATOR
		m_maxRootVictim = 0;
#endif
	}
	assert(!(m_root->previousg));
	if (tmpg.previousg)
		tmpg.previousg->nextg = tmpg.nextg;
	assert(tmpg.nextg);
	tmpg.nextg->previousg = tmpg.previousg;
	if (group.nextg) {
		group.nextg->previousg = &tmpg;
	}
	tmpg.nextg = group.nextg;
	tmpg.previousg = &group;
	tmpg.count = goalVal;
	group.nextg = &tmpg;
	assert(m_root->nextg != m_root);
}

inline void FAST::newGroup(FASTITEM &il, FASTGROUP &group, unsigned int goalVal) {
	assert(group.count < goalVal);
	assert(il.parentg->items != &il);
	FASTGROUP &newgroup = *(m_freegroups[m_gpt++]); //get new group
	newgroup.count = goalVal; // set count to the actual value of the new group
	newgroup.items = &il;
	if (group.nextg) { // if there is another group
		group.nextg->previousg = &newgroup;
	}
	newgroup.nextg = group.nextg;
	newgroup.previousg = &group;
	group.nextg = &newgroup;
	il.parentg = &newgroup;
	assert(m_root->nextg != m_root);
}

int const FAST::size() {
	// This function is no longer correct
	return 0;
	//return sizeof(FAST) + (m_tblsz) * sizeof(FASTITEM*) +
		//(m_nrCounters)*(sizeof(FASTITEM) + sizeof(FASTGROUP) + sizeof(FASTITEM*));
}

inline void FAST::recycleGroup(FASTGROUP & oldgroup)
{
	// Given a group with no counters (oldgroup), we remove it from the SOS and add it back to the pool
	if (oldgroup.nextg) // there is another group
		oldgroup.nextg->previousg = oldgroup.previousg;
	if (oldgroup.previousg)
		oldgroup.previousg->nextg = oldgroup.nextg;
	if (m_root == &oldgroup) // this is the first group
	{
		assert(!(oldgroup.nextg->previousg));
		m_root = oldgroup.nextg;

	}

	//recycle oldgroup.
	m_freegroups[--m_gpt] = &oldgroup;
	// if we have created an empty group, remove it 
	assert(m_root->nextg != m_root);
}

void FAST::clear()
{
	int i;
	for (i = 0; i<m_tblsz; i++)
		m_hashtable[i] = NULL;

	m_groups->count = 0;
	m_groups->nextg = NULL;
	m_groups->previousg = NULL;

	m_groups->items = m_items;// TODO: check if there is more counters to reset
}
