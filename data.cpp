/*
 *  data.cpp
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "data.h"
#include <boost/lexical_cast.hpp> 
#include <fstream>
#include <sstream>
#include <cmath>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>


double Data::mylog(double p1, double base){
	if(p1<=0)
		return 0.0;
	else {
		gsl_sf_result logp1,logbase;
		int status_p1=gsl_sf_log_e(p1,&logp1);
		if(status_p1 != GSL_SUCCESS){
			cout << "Kernel::mylog  GSL Error "<<  status_p1 << " for p1 =" << p1 << endl;
			exit(-1);
		}
		int status_base=gsl_sf_log_e(base,&logbase);  // this division will normalize the entropy value between 0 and 1.
		if(status_base != GSL_SUCCESS){
			cout << "Kernel::mylog GSL Error "<< status_base<< " for base =" << base << endl;
			exit(-1);
		}
		return logp1.val/logbase.val;
	}
}

Data::Data(){
	_mass = 70;
	_rbase=1.263;
	_dt=0.1;
	_N0=400.0;
	_tau=365.0;	
	_frac_csc=0.0;
	_numlsc = 0;
	//are the same accross mamals
	_epsh=0.8476;
	_epsc=0.72;
	_epsi=_epsc;
	_ntimes=91250;
	_age=25;
	_ncompartments=32;
	_stop=10.39;//DP RAT
	_reduction = 3;
	_limit=7;
	_additional=0;
	_threshold = 0.2;
	//for storing data
	_step=100;
	_ofcompartment="compartment.txt";
	_offinal="mammal.txt";
	_ofname = "result.txt";
	_hdlocation ="./";

	//not relevant for the moment
	_epsb=0.90;
	_epsi=0.8476; // what value?
	_p_csc=0;
	_p_imm=0;
	_perc_bound=0.0;
	_treatment=3.0;
	_rcancer = 1.0;
}


void Data::calcFromMass(double mass, double N,double B, double T, double L, double nummonth){ 
	_mass = mass;
	
	double hsc=N * pow(mass, 0.75);
	double temp;
	modf(hsc,&temp); 	
	double diff=hsc-temp;
		
	if(diff >= 0.5)
		_N0 = (temp+1); // always round up to next complete HSC
	else _N0=temp;
	_frac_csc=1/_N0; //always start with 1 LSC
	_numlsc = 1;
	
	_tau = 365.0/(B*pow(mass,-0.25));	
	_dt=T*pow(mass,0.25);	
	
	_age=(L*pow(mass,0.25));
	_ntimes = int((365.0/_dt)*(L*pow(mass,0.25)));
	
	_step=int(nummonth*(30./_dt)); //output per x months
	
	//are the same accross mamals or changd by command line
	_rbase=1.263;
	_epsh=0.8476;
	_epsc=0.72;
	_ncompartments=32;
	_stop=12;//log value at diagnosis
	_reduction = 3; //required log reduction
	_limit=7;
	_additional=0;
	_threshold = 0.2;
	
	//for storing data
	_ofcompartment="iterresult.txt";
	_offinal="mammal.txt";
	_ofname = "result.txt";
	_hdlocation ="./";

	//not relevant for the moment
	_epsb=0.89; //imatinib;
	_epsi=_epsc;
	_p_csc=0;
	_p_imm=0;
	_perc_bound=0.05;
	_treatment=3.0;
	_rcancer = 1.0;
}

Data::Data(const Data& other){
	_dt=other.dt();
	_epsh=other.epsh();
	_epsc=other.epsc();
	_epsb=other.epsb();
	_epsi=other.epsi();
	_p_csc=other.p_csc();
	_p_imm=other.p_imm();
	_frac_csc=other.frac_csc();
	_numlsc = other.numlsc();
	_perc_bound=other.perc_bound();
	_ntimes=other.ntimes();
	_age=other.age();
	_rbase=other.rbase();
	_ncompartments=other.ncompartments();
	_stop=other.stop();
	_reduction=other.reduction();
	_N0=other.N0();
	_tau=other.tau();
	_limit=other.limit();
	_threshold = other.threshold();
	_additional=other.additional();
	_step=other.step();
	_treatment=other.treatment();
	_ofcompartment=other.ofcompartment();
	_offinal=other.offinal();
	_ofname = other.ofname();
	_hdlocation =other.storage();
	_rcancer = other.rcancer();
	_mass = other.mass();
}


ostream& Data::display(ostream& os){
	os << "#Input data for hematopoietic model" << endl;
	os << "\tcalculated :: " << endl;
	os << "\t\tN0 " << _N0 << endl;
	os << "\t\ttau " << _tau << endl;
	os << "\t\tdt " << _dt << endl;
	os << "\t\tntimes " << _ntimes << " iterations = " << _age << " years" << endl;
	os << "\tcommand line :: " << endl;
	os << "\t\tmass " << _mass << endl;
	os << "\t\tlimit " << _limit << endl;
	os << "\t\tstep " << _step << endl;
	os << "\t\tmammal " << _offinal << endl;
	os << "\t\tthreshold " << _threshold << endl;
	os << "\t\tfrac_csc " << _frac_csc << endl;
	os << "\t\tnumlsc " << _numlsc << endl;
	os << "\tfixed :: " << endl;
	os << "\t\trbase " << _rbase << endl;	
	os << "\t\tepsh " << _epsh << endl;
	os << "\t\tespc " << _epsc << endl;
	os << "\t\tespb " << _epsb << endl;
	os << "\t\tespi " << _epsi << endl;
	os << "\t\tp_csc " << _p_csc << endl;
	os << "\t\tp_imm " << _p_imm << endl;
	os << "\t\tperc_bound " << _perc_bound << endl;
	os << "\t\tncompartment " << _ncompartments << endl;
	os << "\t\tadditional " << _additional << endl;
	os << "\t\ttreatment " << _treatment << endl;
	os << "\t\tstop " << _stop << endl;
	os << "\t\treduction " << _reduction << endl;
	os << "\t\trcancer " << _rcancer << endl;
	os << "\tdata collection :: " << endl;	
	os << "\t\toutput 1  " << _ofcompartment << endl;
	os << "\t\toutput 3  " << _ofname << endl;
	os << "\t\tstorage  " << _hdlocation << endl;
	
	return os;
}

