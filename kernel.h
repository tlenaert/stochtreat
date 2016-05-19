/*
 *  kernel.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/5/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __KERNEL_H
#define __KERNEL_H

#include <fstream>
#include "model.h"
#include "reactions.h"
#include "dependency.h"
#include "indexedqueue.h"
#include "rangen.h"

class Kernel {
public:
	Kernel(RanGen& ran, Data& data, unsigned size, double time=0.0):_time(time),_data(data), _dt(data.dt()), _pool(data,size), 
	_depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),_endtime(data.ntimes()*(data.dt()/365.0)),_lsctime(-1){};
	void printAll();
	double execute(RanGen& ran, double t, bool treat);
	bool reachedDiagnosis();
	float getDiagnosis();
	bool reachedReduction();
	float getReduction();
	float whenReduction();
	bool hasLSC();
	float burden();
	void addStochCompSizes(double* data);
	ostream& writeModel(ostream&);
	
private:
	void adjustReactions(RanGen& ran,unsigned compartment, unsigned type);
	bool directMethod(RanGen& ran);
	bool nextMethod(RanGen& ran);
	void detUpdate();
	void treatCells(RanGen& ran);
	
	
	double _time;
	Data _data;
	double _dt;
	Model _pool;
	AllReactions _allr;
	DependencyGraph _depend;
	IndexedQueue _queue;
	double _endtime;
	double _lsctime;
	
};

#endif
