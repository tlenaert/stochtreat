/*
 *  data.h
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Copyright 2006 SWITCH. All rights reserved.
 *
 */


#ifndef __DATA_H
#define __DATA_H

#include <string>
#include <istream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>


class Data {
	public:
		Data(); 
		Data(const Data&); 
		~Data(){};
		
		double dt() const {return _dt;}
		void setdt(double v) {_dt=v;}
		double frac_csc() const {return _frac_csc;}
		void setFrac_csc (double v) {
			_frac_csc = v;
			_numlsc = v*N0();
			double temp;
			modf(_numlsc,&temp); 	
			double diff=_numlsc-temp;
			
			if(diff >= 0.5)
				_numlsc = (temp+1); // always round up to next complete lsc
			else _numlsc=temp;			
		}
		double numlsc() const {return _numlsc;}
		void setNumLsc(double v) {
			_numlsc = v;
			_frac_csc = v/N0();
		} 
		double rbase() const {return _rbase;}
		void setRbase (double v) {_rbase = v;}
		double N0() const {return _N0;}
		void setN0(double v) {_N0 = v;}

                /** Returns the average cell cycle time of HSC.*/
		double tau() const {return _tau;}
                /** Sets the average cell cycle time of HSC.*/
		void setTau(double v) {_tau = v;}

		double mass() const {return _mass;}
		void setMass(double v) {_mass = v;}

		double rcancer() const {return _rcancer;}
		void setRcancer(double v) {_rcancer = v;}
		
		double epsh() const {return _epsh;}
		double epsc() const {return _epsc;}
		double epsb() const {return _epsb;}
		double epsi() const {return _epsi;}

                /** return self-renewal probability epsilon for type. */
		double eps(unsigned type) const {
			switch(type){
				case 0:
					return _epsh;
				case 1:
					return _epsc;
				case 2:
					return _epsi;
				case 3:
					return _epsb;
				default :
					return -1.0;
			}
		}
	
		void setEpsh (double v) {_epsh = v;} 
		void setEpsc (double v) {_epsc = v;} 
		void setEpsb (double v) {_epsb = v;} 
		void setEpsi (double v) {_epsi = v;} 
		void setEps(unsigned type, double v)  {
			switch(type){
				case 0:
					_epsh = v;
					break;
				case 1:
					_epsc = v;
					break;
				case 2:
					_epsi = v;
					break;
				case 3:
					_epsb = v;
					break;
			}
		}
	
		double p_csc() const {return _p_csc;}
		void setPcsc (double v) {_p_csc = v;} 
		double p_imm() const {return _p_imm;}
		void setPimm (double v) {_p_imm = v;} 


                /** returns the treatment rate of all cancer cells. */
		double treatment_rate() const {return _treatment_rate;}

                /** sets the percentage of cells that is affected by treatment. */
		void set_treatment_rate(double v) {_treatment_rate = v;}
		
		
                /** returns the total number of compartments in the model. */
		int ncompartments() const {return _ncompartments;}
                /** Sets the total number of compartments in the model. */
		void setNCompartments(int v) {_ncompartments = v;}

                /** Returns the threshold power of 10 when diagnosis is reached:
                 * When 10^(stop) cells are in the compartment -> diagnosis. */
		double stop() const {return _stop;}

                /** Sets the threshold power of 10 when diagnosis is reached:
                 * When 10^(stop) cells are in the compartment -> diagnosis. */
		void setStop(double v) {_stop = v;}

                /** Returns the recuction level that is required to stop treatment.*/
		double reduction() const {return _reduction	;}
		void setReduction(double v) {_reduction = v;}

		double additional() const {return _additional;}
		void setAdditional(double a)  {_additional=a;}

                /** returns the number of stochstic compartments */
		int nstochcomp() const {return _numstochcomps;}

                /** Sets the number of stochastic compartments.*/
		void setLimit(int l) {_numstochcomps = l;}

		double threshold() const {return _threshold;}
		void setThreshold(double v) {_threshold = v;}
		  

		int step() const {return _step;}
		void setStep(int v) { _step = v;}
                std::string ofcompartment() const {return _ofcompartment;}
		void setOfcompartment (std::string name) {_ofcompartment = name;}
                std::string offinal() const {return _offinal;}
		void setOffinal (std::string name) { _offinal = name;}
		void setOfname(std::string name) {_ofname=name;}
                std::string ofname() const {return _ofname;}

                /** returns treatment time in years.*/
		double treatment_dur() const {return _treatment_duration;}

                /** Returns the maximum time a simulation runs in years.*/
                double getTmax_in_years() const {return _tmax;}

                /** Sets maximum simulation time. */
                void setTmax(double v){_tmax=v;}

                /** Sets treatment time in years that will be used 
                 * in the treatment phase of the simulation. */
		void setTreatment(double t) {_treatment_duration=t;}

                std::string storage() const {return _hdlocation;}
		void setStorage(std::string s) {_hdlocation=s;}
		
                /** Calculates patient parameters from given input.
                 * parameters:
                 * mass     - mass of the modeled animal
                 * Bbase    - log base for the average cell cycle time of hematopeotic stem cells
                 * Sbase    - log base for the deterministic timestep of simulation
                 * Lbase    - log base for maximum simulation time */
		void calcFromMass(double,double,double, double, double, double);

		friend std::ostream & operator<<(std::ostream &o, Data& c){return c.display(o);}
					
	private:
                std::ostream& display(std::ostream&);
                std::istream& read_from_file(std::istream&);

		double _dt; //time step relative to days
		double _epsh; //differentation probability normal cell
		double _epsc; //differentation probability cancer cell
		double _epsb; //differentation probability inhibited cell
		double _epsi; //differentation probability immmune cell
		double _p_csc; //probability that a normal cell turns into a cancer cell
		double _p_imm; //probabilty that a cancer cell turns into an immune cell
		double _frac_csc; //fraction of cancer cells in the stem cell compartment
		double _numlsc; //number of LSC
		double _treatment_rate; //percentage of cells bound to imatinib per day
		double _rbase; //basal metabolic rate
                double _tmax; // maximum simulation time in years
		int _ncompartments;  //number of compartmens in the hematopoeitic system
		double _N0; //Numbe of cells in the stem cell compartment
		double _tau; //maturation time reticulocytes
		int _numstochcomps; //index of first deterministic compartment
		double _additional; //additional number of years to continue simulation after X
		double _treatment_duration; //number of years of treatment
		double _stop; //stop value = diagnosis level
		double _reduction; //stop value (CHANGE now it is the required log reduction in bcr-abl transcript level)
		double _mass;  //mammal mass
		double _rcancer; //difference between replication rates of normal and cancer cells.
		double _threshold; //percentage increase in number of cells for diagnosis

		int _step; 
                std::string _hdlocation; 
		std::string _ofcompartment;
		std::string _offinal;
		std::string _ofname;
};

#endif
