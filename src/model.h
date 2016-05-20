/*
 *  Model.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __Model_H
#define __Model_H

#include <string>
#include <istream>
#include <iostream>
#include "data.h"
#include "rangen.h"

using namespace std;

enum celltypes {H=0,C,I,B}; //H=healthy, C=cancer, I=resitant to TKI and B=Bound to TKI.


class Model {
public:
	Model(Data data, unsigned int size);
	Model(const Model& other);
	~Model(){
		delete[] _compartments;
		delete[] _rates;
		delete[] _previous;
	}
	
	unsigned int numStoch() const {return _numstoch;}
	unsigned int numComp() const {return _numcomp;}
	
	void setRate(unsigned int k, double v);
	double getRate(unsigned k) const;
	
	double getN(unsigned int k) const;
	void setH(unsigned int k, double v); 
	double getH(unsigned int k) const; 
	void setC(unsigned int k, double v); 
	double getC(unsigned int k) const; 
	void setI(unsigned int k, double v); 
	double getI(unsigned int k) const; 
	void setB(unsigned int k, double v); 
	double getB(unsigned int k) const; 

	double retrieveN(unsigned int k) const;
	void storeH(unsigned int k, double v); 
	double retrieveH(unsigned int k) const; 
	void storeC(unsigned int k, double v); 
	double retrieveC(unsigned int k) const; 
	void storeI(unsigned int k, double v); 
	double retrieveI(unsigned int k) const; 
	void storeB(unsigned int k, double v); 
	double retrieveB(unsigned int k) const; 
	
	double lastN() const {return getN(_numcomp-1);}
	double lastH() const {return getH(_numcomp-1);}
	double lastC() const {return getC(_numcomp-1);}
	double lastI() const {return getI(_numcomp-1);}
	double lastB() const {return getB(_numcomp-1);}
	
	double diagRes() const {return _diagnosis;}
	void setDiagRes(double v) {_diagnosis = v;}
	double when() const {return _when;}
	void setWhenReduction(double v) {_when = v;}
	
	friend ostream & operator<<(ostream &o, Model& p){return p.display(o);}
	Model& operator=(const Model&);

	double retrieve(unsigned k, unsigned t);
	void store(unsigned k, unsigned t, double v);

	double get(unsigned k, unsigned t);
	void set(unsigned k, unsigned t, double v);
	void incr(unsigned k, unsigned t, double v);
	void decr(unsigned k, unsigned t, double v);
	
	void memorize();
	bool updateDet(unsigned k, Data& data);
	bool diagnosis(Data& data);
	bool reduction(Data& data);
	float getReduction();
	bool containsLSC();
	bool treatDeterministically(unsigned k, double amount);
	bool treatStochastically(unsigned k, double rate, RanGen& ran);
	
	void calcAlpha();
	double getAlpha() const {return _alpha;}
	
	double diseaseBurden();
	
private:
	double myround(double val);
	double mylog(double p1, double base);

	ostream& display(ostream&);

	unsigned int _endstoch;
	double* _compartments;
	double* _previous;
	double* _rates;
	unsigned _numstoch;
	unsigned _numcomp;
	double _alpha;
	double _diagnosis;
	double _when;
};

#endif
