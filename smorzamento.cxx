#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

double media( vector <double> v);
double somma( vector <double> v);
double stDev( vector <double> v);
int index(vector<double> v, double val);
double somma2_sigma2 (const vector<double> v, const vector<double> s);
double somma_sigma2 (const vector<double> v, const vector<double> s);
double sommaP_sigma2 (const vector<double> x, const vector<double> y, const vector<double> s);
double sigma2 (const vector<double> s);
double delta (const vector<double> v, const vector<double> s);

int main(){

vector <string> file_input = {"dati1.txt","dati2.txt","dati3.txt","dati4.txt"};
vector <string> file_output = {"output_dati1.txt","output_dati2.txt","output_dati3.txt","output_dati4.txt"};
vector <string> analisi_output = {"analisi_dati1.txt","analisi_dati2.txt","analisi_dati3.txt","analisi_dati4.txt"};
vector <string> intervalli_output = {"intervalli_dati1.txt","intervalli_dati2.txt","intervalli_dati3.txt","intervalli_dati4.txt"};
vector <double> stdev_max = {0.04427605,0.0473162,0.0456977,0.0458736};
vector <double> stdev_min = {0.04242585,0.0443194,0.0449084,0.0399874};

//inizio analisi per ogni file

for (int running_file = 0; running_file < file_input.size(); running_file++){


ifstream fin(file_input.at(running_file));
	if (!fin){
		cout << "error " << endl;
		return -1;
}
ofstream fout1(file_output.at(running_file));
	if (!fout1){
		cout << "error " << endl;
		return -
		1;
}
ofstream fout(analisi_output.at(running_file));
	if (!fout){
		cout << "error " << endl;
		return -1;
}
ofstream fout2(intervalli_output.at(running_file));
	if (!fout2){
		cout << "error " << endl;
		return -1;
}


vector <double> tempo0, forzante0, pendolo0, ampiezza, fase;
double t,f,p,a,fa;
while (fin >> t >> f >> p >> a >> fa){
		tempo0.push_back(t);
		forzante0.push_back(f);
		pendolo0.push_back(p*2*M_PI); //angolo in radianti
		ampiezza.push_back(a);
		fase.push_back(fa);	
}

vector <double> tempo, pendolo;
double F0 = 0.001;
for (int i = 1; i < forzante0.size(); i++){
	if (abs(forzante0.at(i-1)) < F0 && abs(forzante0.at(i)) < F0){
		tempo.push_back(tempo0.at(i-1));
		pendolo.push_back(pendolo0.at(i-1));
	}
}		

fout1 << "tempo (s) \t angolo pendolo (rad) " << endl;
for (int i = 0; i < tempo.size(); i++){
	fout1 << tempo.at(i) << "\t\t" << pendolo.at(i) << endl;
	}




//calcolo degli zeri e stima del periodo
vector <double> istanti;
for (int i = 1; i < pendolo.size()-1; i++){
	if (pendolo.at(i-1)*pendolo.at(i) <= 0.0 && pendolo.at(i-1) <= pendolo.at(i)){
		double t = (tempo.at(i-1)+tempo.at(i))/2;
		istanti.push_back(t);
	}
}
//cout << istanti.size() << endl;

int i = 1;
while (i < istanti.size()) {
    if (istanti.at(i) - istanti.at(i-1) < 0.5){
        istanti.erase(istanti.begin() + i); //begin punta al primo elemento + i esimo elemento
    	} 
    else {
        i++;  // Vai avanti solo se non hai cancellato
    }
}

while (istanti.size() > 50){ //con 50 funziona per 1 2, 3 gneh alla fine va sotto lo 0, 4 devo mettere con 40
	istanti.pop_back();
	}
//cout << istanti.size() << endl;
//for (int i = 0; i < istanti.size(); i++){
//	cout << istanti.at(i) << endl;
//	}
	

vector <double> periodi; //prendo solo i pari per scorrelare
for (int i = 0; i < istanti.size()-1; i+=2){
	double T_i = istanti.at(i+1)-istanti.at(i);
	if ( T_i > 0.5){
	periodi.push_back(T_i);
	}
}
//cout << periodi.size() << endl;

//bisogna avere una stima di ordine di grandezza almeno di quanto un dato può fluttuare sulle y
double T = media(periodi);
double sigma_T = stDev(periodi)/sqrt(periodi.size());
fout << "Pseudoperiodo Ts = " << T << " ± " << sigma_T << endl; 
double omegaS = 2*M_PI/T;
double err_omegaS = 2*M_PI*sigma_T/(pow(T,2));
fout << "Pulsazione omega S = " << omegaS << " ± " << err_omegaS << endl;

// Calcolo max e min assoluti
double max = pendolo.at(0);
double min = pendolo.at(0);
for (auto c : pendolo) {
	if (c > max) max = c;
	if (c < min) min = c;
}

int idx_max = index(pendolo, max);
if (idx_max == -1) {
	cerr << "Errore: massimo non trovato!\n";
	return -1;
}
double tmax = tempo.at(idx_max);

int idx_min = index(pendolo, min);
if (idx_min == -1) {
	cerr << "Errore: minimo non trovato!\n";
	return -1;
}
double tmin = tempo.at(idx_min);

//cout << tempo.size() <<"\t"<< pendolo.size() << endl;

//estrazione dei massimi
vector <double> va, vb;
vector<double> massimi;
vector<int> indici_massimi;
for (double i = tmax; i < tempo.at(tempo.size()-1); i += 2 * T) { //tempo.back()
	double a = i - T / 4;
	double b = i + T / 4;
	double max_in_range = 0;
	double index_max = -1;
	va.push_back(a);
	vb.push_back(b);
	for (int j = 0; j < tempo.size(); j++){
		if (tempo.at(j) > a && tempo.at(j) < b){
			if (pendolo.at(j) > max_in_range){
			max_in_range = pendolo.at(j);
			index_max = j;
			}
		}
	}
	
	if (index_max != -1) {
		massimi.push_back(max_in_range);
		indici_massimi.push_back(index_max);
	}
}
//estrazione dei minimi

vector <double> vc, vd;
vector<double> minimi;
vector<int> indici_minimi;
for (double i = tmin; i < tempo.at(tempo.size()-1); i += 2 * T){ //tempo.back()
	double c = i - T / 4;
	double d = i + T / 4;
	double min_in_range = 0;
	double index_min = -1;
	vc.push_back(c);
	vd.push_back(d);
	for (int j = 0; j < tempo.size(); j++){
		if (tempo.at(j) > c && tempo.at(j) < d){
			if (pendolo.at(j) < min_in_range){
			min_in_range = pendolo.at(j);
			index_min = j;
			}
		}
	}
	if (index_min != -1) {
		minimi.push_back(min_in_range);
		indici_minimi.push_back(index_min);
	}
} 


//stampa intervalli ab e cd 

int k = va.size();
if (vb.size() < k){
	k = vb.size();
	}
if (vc.size() < k){
	k = vc.size();
	}
if (vd.size() < k){
	k = vd.size();
	}
	
fout2 << "a\t\tb\t\tc\t\td\t\t" << endl;
for (int i = 0; i < k; i++){
	fout2 << va.at(i) << "\t" << "0" << "\t" << vb.at(i) << "\t" << "0" << "\t" << vc.at(i) << "\t" << "0" << "\t" << vd.at(i) << "\t" << "0" << endl;
	}
	
//calcola istanti max e min
vector <double> tempi_massimi;
for (int i = 0; i < indici_massimi.size(); i++){
	tempi_massimi.push_back(tempo.at(indici_massimi.at(i)));
	}
vector <double> tempi_minimi;
for (int i = 0; i < indici_minimi.size(); i++){
	tempi_minimi.push_back(tempo.at(indici_minimi.at(i)));
	}


//cout << "massimi.size(): " << massimi.size() << endl;
//cout << "minimi.size(): " << minimi.size() << endl;
//cout << "tempi_massimi.size(): " << tempi_massimi.size() << endl;
//cout << "tempi_minimi.size(): " << tempi_minimi.size() << endl;


//offset
vector <double> vector_offsets;
if (minimi.size() <= massimi.size()){
	for (int i = 0; i < minimi.size(); i++){
		vector_offsets.push_back((massimi.at(i)- abs(minimi.at(i)))/2);
	}
}
else if (minimi.size() > massimi.size()){
	for (int i = 0; i < massimi.size(); i++){
		vector_offsets.push_back((massimi.at(i)- abs(minimi.at(i)))/2);
	}
}


double Offset = media(vector_offsets);
double err_offset = stDev(vector_offsets)/sqrt(vector_offsets.size());
fout << "l'Offset è = " << Offset << " ± " << err_offset << endl; //errore piccolo 10% quindi è costante nel nostro t

//scala log prima della correzione dell'offset (serve solo per vedere grafici, NO IN ANALISI)
vector <double> nol_max, nol_min;
for (int i = 0; i < massimi.size(); i++){
	nol_max.push_back(log(massimi.at(i)));
	}
for (int i = 0; i < minimi.size(); i++){
	nol_min.push_back(log(abs(minimi.at(i))));
	}


//correzione dell'offset(errore sistematico)
vector <double> max_off, min_off;
for (int i = 0; i < massimi.size(); i++){
	max_off.push_back(massimi.at(i)-Offset);
	}
for (int i = 0; i < minimi.size(); i++){
	min_off.push_back(minimi.at(i)+Offset);
	}

//errore su ln(max)	
vector <double> errore_ln_massimi;
vector <double> errore_ln_minimi;
for (int i = 0; i < max_off.size(); i++){
	errore_ln_massimi.push_back(stdev_max.at(running_file)/abs(max_off.at(i)));
	}
for (int i = 0; i < min_off.size(); i++){
	errore_ln_minimi.push_back(stdev_min.at(running_file)/abs(min_off.at(i)));
	}

//metto in scala logaritmica
vector <double> ln_massimi, ln_minimi;
for (int i = 0; i < max_off.size(); i++){
	ln_massimi.push_back(log(max_off.at(i)));
	}
for (int i = 0; i < min_off.size(); i++){
	ln_minimi.push_back(log(abs(min_off.at(i))));
	}

//taglio i dati brutti prima dell'interpolazione (per dati 1)
if ( running_file == 0){
while (ln_massimi.size() > 28 && ln_minimi.size() > 28 &&
	errore_ln_massimi.size() > 28 && errore_ln_minimi.size() > 28 &&
	tempi_massimi.size() > 28 && tempi_minimi.size() > 28){
		ln_massimi.pop_back();
		ln_minimi.pop_back();
		errore_ln_massimi.pop_back();
		errore_ln_minimi.pop_back();
		tempi_massimi.pop_back();
		tempi_minimi.pop_back();
}}	
cout << ln_massimi.size() << endl;

//interpolazione max y = ln(tetha omegenea 0) - gamma x
// a = ln(tetha omogenea 0)    b = gamma
double Na_max = somma2_sigma2(tempi_massimi,errore_ln_massimi)*somma_sigma2(ln_massimi,errore_ln_massimi) - somma_sigma2(tempi_massimi,errore_ln_massimi)*sommaP_sigma2(tempi_massimi,ln_massimi,errore_ln_massimi);
double ln_tetha_om_max = Na_max/delta(tempi_massimi,errore_ln_massimi);
double sigma_tetha_om_max = sqrt(somma2_sigma2(tempi_massimi,errore_ln_massimi)/delta(tempi_massimi,errore_ln_massimi));
	
double Nb_max = sigma2(errore_ln_massimi)*sommaP_sigma2(tempi_massimi,ln_massimi,errore_ln_massimi) - somma_sigma2(tempi_massimi,errore_ln_massimi)*somma_sigma2(ln_massimi,errore_ln_massimi);
double gamma_max = Nb_max/delta(tempi_massimi,errore_ln_massimi);
double sigma_gamma_max = sqrt(sigma2(errore_ln_massimi)/delta(tempi_massimi,errore_ln_massimi));
	
double num_max = 0;
for (int i = 0; i < tempi_massimi.size(); i++){
	num_max += pow((ln_massimi.at(i) - ln_tetha_om_max - gamma_max*tempi_massimi.at(i)),2);
	}
double denom_max = tempi_massimi.size() - 2;
double sigmaY_post_max = sqrt(num_max/denom_max);		

fout << "retta interpolante dei massimi y = ln(tetha omegenea 0) - gamma*x" << endl;
fout << "ln(tetha omogenea 0) " << ln_tetha_om_max << " ± " << sigma_tetha_om_max << endl;
fout << "gamma " << abs(gamma_max) << " ± " <<  sigma_gamma_max << endl;
fout << "con sigmaYpost = " << sigmaY_post_max << endl << endl;

// interpolazione min y = ln(tetha omogenea 0) - gamma x
// a = ln(tetha omogenea 0)    b = gamma
double Na_min = somma2_sigma2(tempi_minimi, errore_ln_minimi) * somma_sigma2(ln_minimi, errore_ln_minimi) 
                - somma_sigma2(tempi_minimi, errore_ln_minimi) * sommaP_sigma2(tempi_minimi, ln_minimi, errore_ln_minimi);
double ln_tetha_om_min = Na_min / delta(tempi_minimi, errore_ln_minimi);
double sigma_tetha_om_min = sqrt(somma2_sigma2(tempi_minimi, errore_ln_minimi) / delta(tempi_minimi, errore_ln_minimi));

double Nb_min = sigma2(errore_ln_minimi) * sommaP_sigma2(tempi_minimi, ln_minimi, errore_ln_minimi) 
                - somma_sigma2(tempi_minimi, errore_ln_minimi) * somma_sigma2(ln_minimi, errore_ln_minimi);
double gamma_min = Nb_min / delta(tempi_minimi, errore_ln_minimi);
double sigma_gamma_min = sqrt(sigma2(errore_ln_minimi) / delta(tempi_minimi, errore_ln_minimi));

double num_min = 0;
for (int i = 0; i < tempi_minimi.size(); i++) {
    num_min += pow((ln_minimi.at(i) - ln_tetha_om_min - gamma_min * tempi_minimi.at(i)), 2);
}
double denom_min = tempi_minimi.size() - 2;
double sigmaY_post_min = sqrt(num_min / denom_min);

fout << "retta interpolante dei minimi y = ln(tetha omogenea 0) - gamma*x" << endl;
fout << "ln(tetha omogenea 0) " << ln_tetha_om_min << " ± " << sigma_tetha_om_min << endl;
fout << "gamma " << abs(gamma_min) << " ± " << sigma_gamma_min << endl;
fout << "con sigmaYpost = " << sigmaY_post_min << endl << endl;

//stima di gamma
double gamma = abs((gamma_max + gamma_min)/2);
double err_gamma = sqrt( pow(sigma_gamma_max,2) + pow(sigma_gamma_min,2))/2;
fout << "la miglior stima di gamma è = " << gamma << " ± " << err_gamma << endl;

//stima tetha omogenea
double tetha_om = exp((ln_tetha_om_max + ln_tetha_om_min)/2);
double err_tetha = tetha_om*(sqrt( pow(sigma_tetha_om_max,2) + pow(sigma_tetha_om_min,2))/2); 
fout << "la miglior stima di tetha omogenea è = " << tetha_om << " ± " << err_tetha << endl << endl;

//calcolo omegaR attesa
double omegaR = sqrt(omegaS*omegaS - gamma*gamma);
double de_omegaS = omegaS/(sqrt(omegaS*omegaS-gamma*gamma));
double de_gamma = -gamma/(sqrt(omegaS*omegaS-gamma*gamma));
double err_omegaR = sqrt(pow(de_omegaS,2)*pow(err_omegaS,2) + pow(de_gamma,2)*pow(err_gamma,2));
fout <<"la miglior stima di omegaR attesa è = " << omegaR << " ± " << err_omegaR << endl << endl;

//calcolo omegaO
double omegaO = sqrt(omegaS*omegaS + gamma*gamma);
double De_omegaS = omegaS/(sqrt(omegaS*omegaS+gamma*gamma));
double De_gamma = gamma/(sqrt(omegaS*omegaS+gamma*gamma));
double err_omegaO = sqrt(pow(De_omegaS,2)*pow(err_omegaS,2) + pow(De_gamma,2)*pow(err_gamma,2));
fout <<"la miglior stima di omegaO è = " << omegaO << " ± " << err_omegaO << endl << endl;

fout << "t max (s) \t max (rad)\t max corretti \t errore max \t ln(max) \t ln(max) corretti \t errore ln(max)" << 
	"\tt min (s) \tmin (rad) \t min corretti \t errore min \t ln(min) \t ln(min) corretti \t errore ln(max)"<< endl;
	

if (running_file == 0){
	for ( int i = 0; i < ln_massimi.size(); i++){
		fout << tempi_massimi.at(i) << "    " << massimi.at(i) 
		<< "    " << max_off.at(i) << "    " << stdev_max.at(running_file) << "    " << nol_max.at(i)
		<< "    " << ln_massimi.at(i) << "    " << errore_ln_massimi.at(i) 
		<< "    " << tempi_minimi.at(i) << "    "  << minimi.at(i) 
		<< "   " << min_off.at(i) << "   " << stdev_min.at(running_file) << "    " << nol_min.at(i) 
		<< "   " << ln_minimi.at(i) << "    " << errore_ln_minimi.at(i) << endl;
	}
}
else{
	if (minimi.size() > massimi.size()){
	for (int i = 0; i < massimi.size(); i++){
		fout << tempi_massimi.at(i) << "    " << massimi.at(i) 
		<< "    " << max_off.at(i) << "    " << stdev_max.at(running_file) << "    " << nol_max.at(i)
		<< "    " << ln_massimi.at(i) << "    " << errore_ln_massimi.at(i) 
		<< "    " << tempi_minimi.at(i) << "    "  << minimi.at(i) 
		<< "   " << min_off.at(i) << "   " << stdev_min.at(running_file) << "    " << nol_min.at(i) 
		<< "   " << ln_minimi.at(i) << "    " << errore_ln_minimi.at(i) << endl;
		}
	}

	else if (minimi.size() <= massimi.size()){
	for (int i = 0; i < minimi.size(); i++){
		fout << tempi_massimi.at(i) << "    " << massimi.at(i) 
		<< "    " << max_off.at(i) << "    " << stdev_max.at(running_file) << "    " << nol_max.at(i)
		<< "    " << ln_massimi.at(i) << "    " << errore_ln_massimi.at(i) 
		<< "    " << tempi_minimi.at(i) << "    "  << minimi.at(i) 
		<< "   " << min_off.at(i) << "   " << stdev_min.at(running_file) << "    " << nol_min.at(i) 
		<< "   " << ln_minimi.at(i) << "    " << errore_ln_minimi.at(i) << endl;
		}
	}
}


}

return 0;
}


double media (vector<double> v){
	double sum = 0;
	for (auto c : v){
		sum += c;
		}
	double media = sum/v.size();
	return media;
}

double somma(const vector<double> v){
    double sum = 0;
    for (auto c : v) {
        sum += c;
    }
    return sum;
}

double stDev (vector<double> v){
    double m = media(v);
    double s = 0;
    for (auto c : v) {
        s += (c - m) * (c - m);
    }
    return sqrt(s / (v.size() - 1));
}

int index(const vector<double> v, double val){
	for (int i = 0; i < v.size(); i++){
		if (v.at(i) == val){
			return i;
		}
	}
	return -1; // valore non trovato
}

//funzioni per interpolazione

double somma2_sigma2 (const vector<double> v, const vector<double> s){
	double sum = 0;
	for (int i = 0; i < v.size(); i++){
		sum += (v.at(i)*v.at(i))/(s.at(i)*s.at(i));
		}
	return sum;
}

double somma_sigma2 (const vector<double> v, const vector<double> s){
	double sum = 0;
	for (int i = 0; i < v.size(); i++){
		sum += v.at(i)/(s.at(i)*s.at(i));
		}
	return sum;
}	
		
double sommaP_sigma2 (const vector<double> x, const vector<double> y, const vector<double> s){
	double sum = 0;
	for (int i = 0; i < x.size(); i++){
		sum += (x.at(i)*y.at(i))/(s.at(i)*s.at(i));
		}
	return sum;
}

double sigma2 (const vector<double> s){
	double sum = 0;
	for (int i = 0; i < s.size(); i++){
		sum += 1/(s.at(i)*s.at(i));
		}
	return sum;
}	

double delta (const vector<double> v, const vector<double> s){
	double t1 = sigma2(s)*somma2_sigma2(v,s);
	double t2 = somma_sigma2(v,s)*somma_sigma2(v,s);
	double delta = t1 - t2;
	return delta;
}




