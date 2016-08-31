/*
 *  reactions.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/1/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "reactions.h"

Reaction::Reaction(const Reaction& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
}

Reaction& Reaction::operator=(const Reaction& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	return *this;
}

SelfRenewal& SelfRenewal::operator=(const SelfRenewal& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();	return *this;
}


Differentation& Differentation::operator=(const Differentation& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	return *this;
}

StemCellRenewal& StemCellRenewal::operator=(const StemCellRenewal& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	_ran = other.ran(); 
	return *this;
}

MoranReaction::MoranReaction(const MoranReaction& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_reactant2=other.reactant2Comp();
	_reactant2=other.reactant2();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
}


MoranReaction& MoranReaction::operator=(const MoranReaction& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_reactant2=other.reactant2Comp();
	_reactant2=other.reactant2();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	return *this;
}


MoranRenewal& MoranRenewal::operator=(const MoranRenewal& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_reactant2=other.reactant2Comp();
	_reactant2=other.reactant2();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	return *this;
}



MoranDifferentiation& MoranDifferentiation::operator=(const MoranDifferentiation& other){
	_comp_reactant=other.reactantComp();
	_reactant=other.reactant();	
	_comp_reactant2=other.reactant2Comp();
	_reactant2=other.reactant2();	
	_comp_product1=other.product1Comp();	
	_product1=other.product1();	
	_comp_product2=other.product2Comp();	
	_product2=other.product2();	
	_r=other.rate();
	_a=other.propensity();
	_ptime = other.putativeTime();
	_lasta = other.last_a();
	_lastptime = other.last_ptime();
	_time_zero = other.getTZero();
	_intype=other.inType();
	_ingraph=other.inGraph();
	_used =other.used();
	return *this;
}


ostream& Reaction::display(ostream& os){
	os << "{";
	os <<_comp_reactant<<":" <<_reactant<<" -> ";
	os <<_comp_product1<<":" <<_product1<<" + ";
	os <<_comp_product2<<":" <<_product2;
	os << "}[";
	os << _a <<" - " << _intype << " - " << _ingraph;
	os << "]"<< "(" << _used << ")";
	return os;
}

ostream& MoranReaction::display(ostream& os){
	os << "{";
	os <<_comp_reactant<<":" <<_reactant<<" + ";
	os <<_comp_reactant2<<":" <<_reactant2<<" -> ";
	os <<"2*"<<_comp_product1<<":" <<_product1<<" + ";
	os <<"2*"<<_comp_product2<<":" <<_product2;
	os << "}[";
	os << _a <<" - " << _intype << " - " << _ingraph;
	os << "]" << "(" << _used << ")";
	return os;
}

bool Reaction::apply(Model& pool, double time){
	//Reaction format : Tp(k) -> Tp(k') + Tp(k'')
	//decrease reactants, increase products
//	cout << "before[" << reactantComp() << ", " << reactant() << ", " << pool.get(reactantComp(), reactant()) << "]" << endl;
	if(pool.retrieve(reactantComp(), reactant()) > 0){
//		if(reactant() == C) cout << "#reaction on C " << pool.getC(reactantComp()) << endl; 
		pool.decr(reactantComp(), reactant(),1);
		pool.incr(product1Comp(), product1(),1);
		pool.incr(product2Comp(), product2(),1);
//		cout << "after[" << reactantComp() << ", " << reactant() << ", " << pool.get(reactantComp(), reactant()) << "]" << endl;
		this->incrUsed();
		return false; 
	}
	return false;
}

bool StemCellRenewal::apply(Model& pool, double time){
	if(pool.retrieve(0, reactant()) > 0){
		bool lsc_moved = false;
		pool.decr(0, reactant(),1);
		pool.incr(0, product1(),1);
		pool.incr(0, product2(),1);
		//but now pick randomly a cell in productComp and remove ti to the next compartment
		double tempC = pool.getC(0);
		double tempH = pool.getH(0);
		double tempI = pool.getI(0);
		double tempN = pool.getN(0);
		double elm=_ran->ranval(0,int(tempN));
		
		if (tempC > 0 && elm < tempC) {
			pool.decr(0, C, 1);
			pool.incr(1, C, 1);
			lsc_moved=true;
		}
		else if(tempH > 0 && elm >=tempC && elm < (tempC+tempH)){ 
			pool.decr(0, H, 1);
			pool.incr(1, H, 1);
		}
		else if(tempI > 0 && elm >=(tempC+tempH) && elm < (tempC+tempH+tempI)){ 
			pool.decr(0, I, 1);
			pool.incr(1, I, 1);
		}		
		this->incrUsed();
		return lsc_moved; 
	}
	return false;
}

//bool StemCellRenewal::apply(Model& pool, double time){
//	if(pool.retrieve(0, reactant()) > 0 && reactant() == H){
//		bool lsc_moved = false;
//		pool.decr(0, reactant(), 1);
//		pool.incr(0, product1(), 1);
//		pool.incr(0, product2(), 1);
//		
//		if (pool.getC(0) > 0) {
//			pool.decr(0, C, 1);
//			pool.incr(1, C, 1);
//			lsc_moved=true;
//		}
//		else {
//			pool.decr(0, H, 1);
//			pool.incr(1, H, 1);
//		}
//
//		this->incrUsed();
//		return lsc_moved; 
//	}
//	return false;
//}


bool MoranReaction::apply(Model& pool, double time){
	//Reaction format : Tp(k) + oTp(k) -> 2Tp(k') + 2oTp(k'')
	//decrease reactants, increase products
//	cout << "before1[" << reactantComp() << ", " << reactant() << ", " << pool.get(reactantComp(), reactant()) << "]\t";
//	cout << "before2[" << reactant2Comp() << ", " << reactant2() << ", " << pool.get(reactant2Comp(), reactant2()) << "]" << endl;

	if( (pool.retrieve(reactantComp(), reactant()) > 0) && (pool.retrieve(reactant2Comp(), reactant2()) > 0) ){
		pool.decr(reactantComp(), reactant(),1);
		pool.decr(reactant2Comp(), reactant2(),1);
//		cout << "after1[" << reactantComp() << ", " << reactant() << ", " << pool.get(reactantComp(), reactant()) << "]\t";
//		cout << "after2[" << reactant2Comp() << ", " << reactant2() << ", " << pool.get(reactant2Comp(), reactant2()) << "]" << endl;

		pool.incr(product1Comp(), product1(),2);
		pool.incr(product2Comp(), product2(),2);
//		cout << "after1[" << product1Comp() << ", " << product1() << ", " << pool.get(product1Comp(), product1()) << "]\t";
//		cout << "after2[" << product2Comp() << ", " << product2() << ", " << pool.get(product2Comp(), product2()) << "]" << endl;
		return false; 
	}
	return false;
}

bool Reaction::sufficientReactants(Model& pool){
	return (pool.get(reactantComp(), reactant()) > 0);
}

double Reaction::reactantFactor(Model& pool){
	return (double)pool.get(reactantComp(), reactant());
}

bool MoranReaction::sufficientReactants(Model& pool){
	if(reactant() != reactant2()){
		return  ( (pool.get(reactantComp(), reactant()) > 0) && (pool.get(reactant2Comp(), reactant2()) > 0) );
	}
	return (pool.get(reactantComp(), reactant()) >= 1);
}

//double MoranReaction::reactantFactor(Model& pool){
//	return (double)(pool.get(reactantComp(), reactant()));
//}

unsigned int AllReactions::add(Reaction* r){
	int loc = (int)_all.size();
	_all.push_back(r);
	return loc;
}

Reaction* AllReactions::operator[](unsigned pos){
	return _all[pos];
}
