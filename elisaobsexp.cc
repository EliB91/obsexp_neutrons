#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>

// Toolkit inclundes
#include "STClibrary.h"
#include "../events.h"
#include "../catalogs.h"
#include "../maptools.h"
#include "../projmap.h"
#include "../agntools.h"
#include "../scanner.h"
#include "../coverage.h"
#include "../common.h"

// ROOT
#include "TRint.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TStyle.h"
using namespace std;

int main(){

////////////////////////////////////////////////////////////////
//                                                            //
//              EVENTI OSSERVATI                              //
//                                                            //
////////////////////////////////////////////////////////////////

// PRIMO PASSO: SCEGLIERE IL NUMERO DI EVENTI NEL FILE //

//1)PRIMO STUDIO DATI CON FILE ENERGIE PAPER APJ2014
/*int h12 = 621375;
int h23 = 135444;
int h3 = 97451;
int h1 = 854270; //eventi_tot
 INFILL - NON ANCORA FATTI
int i1 = 1808183;
int i2 = 636106;
int i3 = 172287;
int i4 = 34534;
int itot = 2651110;*/

    //2)Nuova analisi dati con file Energie Neutrons1EeV
    int h1 = 1847823; //eventi_totE
    int h12 = 1439459;
    int h23 = 255048;
    int h3 = 153316;



////////////////////////////////
int eventi = h3;
////////////////////////////////


ifstream inn;

// CAMBIARE CATALOGO IN BASE ALLE SORGENTI E ALL'ENERGIA //
inn.open("/Users/Brisa/Documenti/cat/E_catalogfin/2f_3.txt"); //devo mettere i final calcolati con la risoluzione angolare calcolata da Danelise
// CAMBIARE CATALOGO IN BASE ALLE SORGENTI E ALL'ENERGIA //

int dim;
inn >> dim;

double * ra_target = new double[dim];
double * decl_target = new double[dim];
double * angres = new double[dim];

int * obs_target = new int[dim];

for(int a=0; a<dim; a++){
    inn >> ra_target[a];
    inn >> decl_target[a];
    inn >> angres[a];
};
inn.close();

ifstream in;


    //QUA NEL SECONDO STEP STO RIFACENDO I CALCOLI USANDO IL NUOVO FILE DI EVENTI NeutronsnewE

// CAMBIARE FILE DI EVENTI IN BASE A ENERGIA E INFILL //

in.open("/Users/Brisa/Documenti/NeutronsnewE/eventi_3E.txt");
// CAMBIARE FILE DI EVENTI IN BASE A ENERGIA E INFILL //

// LE COLONNE SONO LE STESSE: Auger Id - Local Theta - dTheta - Local Phi - dPhi - L -B -Ra - Dec - UTC Time - Tcore - Energy

double * ra_event = new double[eventi];
double * decl_event = new double[eventi];
double * en_event = new double[eventi];
double e; // variabile appoggio per colonne

for(int i = 0; i<eventi; i++){
    in >> e;
    in >> e;
    in >> e;
    in >> e;
    in >> e;
    in >> e;
    in >> e;
    in >> ra_event[i];
    in >> decl_event[i];
    in >> e;
    in >> e;
    in >> en_event[i];
};

in.close();

ofstream out;



// CAMBIARE NOME FILE DI OUTPUT //
out.open("/Users/Brisa/Documenti/analisi/obs/E_obs_cat2f_3.cat");
// CAMBIARE NOME FILE DI OUTPUT //



for(int p = 0; p<dim; p++){

    double conta=0;
    //decide se un evento è o no nel target: la ra di una sorgente candidata in coordinate equat è (a,d) e il suo raggio target è 1,05*ra //
    for(int y = 0; y<eventi; y++){
        double adist = AngularDistance(ra_target[p], decl_target[p], ra_event[y], decl_event[y]);
        if(adist < angres[p]){
            obs_target[p] = obs_target[p] +1;
        }

    }
    cout << ra_target[p] << "         "<< decl_target[p] << "           " << obs_target[p]  << endl;
    out << ra_target[p] << "     "<< decl_target[p] << "           " << obs_target[p]  << endl;
    }

    out.close();


    ////////////////////////////////////////////////////////////////
    //                                                            //
    //              SIMULAZIONE EVENTI ATTESI                     //
    //                                                            //
    ////////////////////////////////////////////////////////////////


	int * cont = new int[dim] ; //vettore per contare il numero di eventi simulati riconducibili ad ogni sorgente in ogni singolo data set
	int * total = new int[dim]; // vettore per contare il numero totale di eventi simulati riconducibili ad ogni sorgente nel totale dei data set

	for(int b=0; b<dim; b++){ //inizializzo a zero i vettori di appoggio per contare gli eventi
		cont[b]=0;
		total[b]=0;
	};



	vector<TEvent> events;


    // CAMBIARE NOME DEL FILE //
	string eventFile = "/Users/Brisa/Documenti/NeutronsnewE/eventi_3E.txt"; //per quelli 1 devi scrivere eventi_totE
     // CAMBIARE NOME DEL FILE //


	cout << "Reading events file " << eventFile << endl;
      	events = GetEvents(eventFile);
      	long nEvents = events.size();
      	if( !events.size() ) {cout << "Program Failed: No events read. Exiting." << endl;
	 exit(0);
	} else {
	cout << "input file OK : " << nEvents << " events read." << endl << endl;
	}

	double latSite = kConstantsTK::AugerSouthLatitude;
  	double lonSite = kConstantsTK::AugerSouthLongitude;
  	double thetaMax = 60.;
	string binningType = "EVENTS";
  	string scramblingType = "UTC+JD";
    unsigned int nBins = 6;

    double dist;

    //PER LO SCRAMBLING NON CAMBIA MOLTO - PER Etot 1000
    //per gli altri 100

    // CAMBIARE NUMERO DI INTERAZIONI //
    int iter = 100;
    // CAMBIARE NUMERO DI INTERAZIONI //


 //fa lo scrambling dei dati, ridistribuisce in modo casuale del tempo di arrivo di ciascun evento mantenendo E,zenith, azimuth
for(int y=0; y<iter; y++){ //ciclo la simulazione per un numero n di data set
	vector<TEvent> EventiRandom = ScrambleData(events, nBins, binningType, latSite, lonSite, thetaMax, scramblingType);

	for(int i = 0; i<dim; i++){ //ciclo su ogni sorgente
		dist = 0; //inizializzo nulla la distanza evento-target (superfluo)
		cont[i] = 0; // in ogni data set simulato, per ogni sorgente inizializzo a zero il contatore degli eventi associati
		for(int j=0; j<nEvents; j++){ //ciclo su ogni evento del y-esimo data set
			if(isnan(EventiRandom[j].fRa)||  isnan(EventiRandom[j].fDec)){ //if di controllo se un evento è stato generato con errori //comando isnan=is Not a number
					j--;
					 continue;
			}
			dist = AngularDistance(ra_target[i], decl_target[i], EventiRandom[j].fRa, EventiRandom[j].fDec);

			if(dist < angres[i]){

				cont[i]++;

			}
		}

	total[i] = total[i] + cont[i]; //conto gli eventi totali (serve per fare la media) sommando per ogni data set, per ogni sorgente, il numero di eventi associati, total, differentemente da cont non si azzera mai così tiene per ogni sorgente il conteggio totale

	}
cout << y+1 << endl;
};

	ofstream out2;




    // CAMBIARE NOME DI FILE DI OUTPUT //
	out2.open("/Users/Brisa/Documenti/analisi/exp/E_exp_cat2f_3_100.cat");
    // CAMBIARE NOME DI FILE DI OUTPUT //




	double *exp = new double[dim];
	cout << "R.A.[°]          Decl.[°]            Exp(total)           Exp "<< endl;
    for(int c = 0; c<dim; c++){


        //ATTENZIONE CAMBIARE N INTERAZIONE//


	exp[c] = total[c]/100.; //effetuo la media tra il numero degli eventi attesi in ogni data set simulato


	cout << ra_target[c]<<"         "<<decl_target[c]<<"                    " <<total[c] <<  "              " << exp[c] <<  endl;
	out2 << ra_target[c]<<"         "<<decl_target[c]<<"                    " <<total[c] <<  "              " << exp[c] << endl;

	}

    delete [] cont;
    delete [] total;
    delete [] exp;
    delete [] ra_event;
    delete [] decl_event;
    delete [] en_event;
    delete [] ra_target;
    delete [] decl_target;
    delete [] angres;
    delete [] obs_target;

return 0;
}

