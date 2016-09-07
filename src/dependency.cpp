/*
 *  dependency.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "dependency.h"

DependencyNode::DependencyNode(const DependencyNode& other){
	_affecting.clear();
	_idx = other.reaction();
	vector<DependencyNode*>::const_iterator start = other.begin();
	vector<DependencyNode*>::const_iterator stop = other.end();
	while(start !=stop){
		affects(*start);
		++start;
	}
}

DependencyNode& DependencyNode::operator=(const DependencyNode& other){
	_affecting.clear();
	_idx =other.reaction();
	vector<DependencyNode*>::const_iterator start = other.begin();
	vector<DependencyNode*>::const_iterator stop = other.end();
	while(start !=stop){
		affects(*start);
		++start;
	}
	return *this;
}

ostream& DependencyNode::display(ostream& os){
	os << "<" << _idx << "|"<< _affecting.size()<<"|";
	for(unsigned i=0; i< _affecting.size(); ++i){
		unsigned tmp =  _affecting[i]->reaction();
		os << tmp;
		if(i < (_affecting.size() -1))
			os << ",";
	}
	os << ">";
	return os;
}



//need to add a mutation reaction at a later stage; mutations occur between cell types. 
DependencyGraph::DependencyGraph(Model& pool, Data& data, AllReactions& all, RanGen& ran){
	//for every type of cell in every compartment there is one self-renewal and one differentiation reaction
	//note that compartment 0 is a special case, since it's size needs to remain constant
	unsigned pos;
	double sum=0.0;
	unsigned numtypes = 4;
	for(unsigned k=1; k < pool.numStoch(); ++k) { //every compartment
		for(unsigned tp = 0; tp < numtypes; ++tp) {
//			cout << "(1-eps) = " << (1.0-data.eps(tp)) << "\t rate= " << pool.getRate(k) << endl; 
			SelfRenewal *self= new SelfRenewal(k,tp,(1.0 - data.eps(tp))*pool.getRate(k));// TP(k) -> TP(k)+TP(k)
			self->setPropensity((self->sufficientReactants(pool)?self->reactantFactor(pool):0.0)); 
			sum += self->propensity();
			pos=all.add(self);
			DependencyNode *selfnode=new DependencyNode(pos);
			
			Differentation *diff= new Differentation(k,tp,data.eps(tp)*pool.getRate(k));// TP(k) -> TP(k+1)+TP(k+1)
			diff->setPropensity((diff->sufficientReactants(pool)?diff->reactantFactor(pool):0.0)); 
			sum += diff->propensity();
			pos=all.add(diff);
			DependencyNode *diffnode=new DependencyNode(pos);
			
			//add dependencies between both nodes and themselves
			selfnode->affects(selfnode);
			selfnode->affects(diffnode);
			diffnode->affects(diffnode);
			diffnode->affects(selfnode);
			
			//both nodes are also affected by certain reactions from the previous compartment.			
			//need to add them. Carefull with the reactions in compartment 0.
			if( (k-1) > 0){ //starting from k=2
				_diffnodes[((k-2)*4 + tp)]->affects(selfnode);
				_diffnodes[((k-2)*4 + tp)]->affects(diffnode);
			}
			//add both nodes to their corresponding vectors.
			pos=add(SELF_RENEWAL,selfnode);
			self->setDG(SELF_RENEWAL,pos);
			pos=add(DIFFERENTATION,diffnode);
			diff->setDG(DIFFERENTATION, pos);	
//			cout << "Self renewal, K=" << k << ", reaction " << *self << endl;
//			cout << "Differentiation, K=" << k << ", reaction " << *diff << endl;
		}
	}
	//add reactions for the stem cell compartment:
	for(int first_tp = 0; first_tp < 3; ++first_tp) { // excluding the bound type, i.e.tp = 3
		StemCellRenewal *scr=new StemCellRenewal(&ran,first_tp,pool.getRate(0));
		scr->setPropensity(scr->sufficientReactants(pool)?scr->reactantFactor(pool):0.0);
		sum += scr->propensity();
		pos=all.add(scr);
		DependencyNode *scrnode=new DependencyNode(pos);
		pos=add(MORAN,scrnode);
		scr->setDG(MORAN, pos);
//		cout << "# stem cell renewal : " << *scr << "\t propensity " << <<endl;

		scrnode->affects(scrnode);
//		cout << "#added loop dependency \n";
		
		if (_selfnodes.size() > 0 && _diffnodes.size() > 0){
			for(int second = 0; second < 3; ++second) { 
				DependencyNode *nextself, *nextdiff;
				nextself = get(SELF_RENEWAL,second); 
				nextdiff = get(DIFFERENTATION,second);
				scrnode->affects(nextself);
//				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextself->reaction()]) << endl;
				scrnode->affects(nextdiff);
//				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[nextdiff->reaction()]) << endl;
			}
//			cout << "#added all dependencies between stem cell reactions and those in the next compartment\n";
		}

	}

	for(unsigned first_tp = 0; first_tp < 3; ++first_tp) { 
		DependencyNode *scrnode =  get(MORAN,first_tp);
		for(unsigned second_tp = 0; second_tp < 3; ++second_tp) {
			if(first_tp != second_tp){
				DependencyNode *oscrnode =  get(MORAN,second_tp);
				scrnode->affects(oscrnode);
//				cout << *(all[scrnode->reaction()]) << "is linked to " << *(all[oscrnode->reaction()]) << endl;

			}
		}
	}
	all.setPropSum(sum);
//	cout << "#finished creation DependencyGraph, number of reactions = "<< all.size() << endl;
}


ostream& DependencyGraph::display(ostream& os){
	os << "DependencyGraph [" << endl << "Moran reactions";
	vector<DependencyNode*>::iterator startmoran = _morannodes.begin();
	vector<DependencyNode*>::iterator stopmoran = _morannodes.end();
	while(startmoran !=stopmoran){
		os << endl << **startmoran;
		++startmoran;
	}
	os << endl<< "Selfrenewal";
	vector<DependencyNode*>::iterator startself = _selfnodes.begin();
	vector<DependencyNode*>::iterator stopself = _selfnodes.end();
	while(startself !=stopself){
		os << endl << **startself;
		++startself;
	}
	os << endl << "Differentiation";
	vector<DependencyNode*>::iterator startdiff = _diffnodes.begin();
	vector<DependencyNode*>::iterator stopdiff = _diffnodes.end();
	while(startdiff !=stopdiff){
		os << endl << **startdiff;
		++startdiff;
	}
	os << "]";
	return os;
}




