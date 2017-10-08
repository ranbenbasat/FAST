/********************************************************************
Implementation of the Fast Flow Volume Estimation paper (ACM ICDCN 2018) 
by R. Ben Basat, G.Einziger and R.Friedman
Implementation by R.Ben Basat and G.Einziger
Last modified: 2017-10
*********************************************************************/
#pragma once
#ifndef FAST_h
#define FAST_h
#include "lossycount.h" // Stealing some defs from Cormode's
#define FAST_HASHMULT 3
#define UNDER_ESTIMATOR 0

typedef struct FAST_item FASTITEM;
typedef struct FAST_group FASTGROUP;

//typedef long long LCUWT;
struct FAST_group
{
	LCUWT count;
	FASTITEM *items;
	FASTGROUP *previousg, *nextg;
};

struct FAST_item
{
	unsigned int item;
	int hash;
	LCUWT delta;
	int remainder;
	FASTGROUP *parentg;
	FASTITEM *previousi, *nexti;
	FASTITEM *nexting, *previousing;
};

class FAST {
public:
	FAST(float fPhi, int M, float gamma); //C'tor - gets the error parameter, max size and gamma
	~FAST();
	void update(int item, int weight); // Processes a packet of flow (item) with size (weight)
	int const size();
	void clear();
	unsigned int const query(unsigned int item); // Estimates the size of flow (item)
private:
	void recycleGroup(FASTGROUP& oldgroup);
	inline void connectToGroup(FASTITEM &newi, FASTGROUP &tmpg);
	void const showGroups();
	void insertIntoHashtable(FASTITEM &newi, int i, unsigned int newitem);
	FASTITEM & takoverMinimal(unsigned int item, int h);
	FASTITEM * getCounter(unsigned int item, int h);
	FASTITEM * getNewCounter(unsigned int nrIncs);
	inline int const getHash(unsigned int item);
	void advanceCounter(FASTITEM &il, unsigned int nrIncs);
	void removeFromHash(FASTITEM & il);
	FASTGROUP & getLastGroup(FASTGROUP &oldgroup, unsigned int goalVal);
	void putInNewGroup(FASTITEM &newi, FASTGROUP & tmpg);
	void moveGroup(FASTITEM &il, FASTGROUP &group, unsigned int goalVal);
	void addToGroup(FASTGROUP & group, FASTITEM & il);
	void newGroup(FASTITEM &il, FASTGROUP &group, unsigned int goalVal);
	void disconnectCounter(FASTITEM & il);
	void computeNrIncsAndAdvance(int newVal, FASTITEM &il);
	//LCUWT m_n;
	int m_gpt;
	int m_tblsz;
	int m_maxRootVictim;
	long long m_a, m_b;
	FASTGROUP * m_root;
	int m_stepSize;
	int m_stepSizeMinusOne;
	int m_M;
	float m_gamma;
	int m_nrCounters;
	unsigned int m_nrFree;
#ifdef LCU_SIZE
	LCUITEM items[LCU_SIZE];
	LCUGROUP groups[LCU_SIZE];
	LCUGROUP *freegroups[LCU_SIZE];
	LCUITEM* hashtable[LCU_TBLSIZE];
#else
	FASTITEM  *m_items;
	FASTGROUP *m_groups;
	FASTGROUP **m_freegroups;
	FASTITEM  **m_hashtable;

#endif
};

#endif
