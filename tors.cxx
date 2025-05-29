#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib> // Per system()
#include <cstdio>  // Per popen(), pclose()

using namespace std;

double media(vector<double> v);
double stDev(vector<double> v);
double somma(const vector<double> v);
double min(vector<double> v);
double max(vector<double> v);
vector<pair<double, double>> dimezza(const vector<pair<double, double>>& v);
double mediap (const vector<pair<double, double>>& v);
double stDevp (const vector<pair<double, double>>& v);
\\ credo che dovro fare altre funzioni per la miglior stima di wr
struct risonanza {
    vector <double> tempo ;
    vector <double> forzante ;
    vector <double> pendolo ;
    vector <double> ampiezza ;
    vector <double> fase ;
};


int main(){

vector <double> ForzForzante = 
{ 0.96 , 0.959 , 0.971 , 0.972 ,
 0.6955 , 0.9419 , 0.9476 , 0.9533 ,
 0.9562 , 0.9614 , 0.9648 , 0.9655 ,
 0.9696 , 0.9751 , 0.9781 , 0.9841 ,
 0.9901 , 0.93619 , 0.96275 , 0.96413 , 
 0.96619 , 0.96688 , 0.96825 , 0.96825 };

vector <string> file_input =
{"FF096.txt","FF0959.txt","FF0971.txt","FF0972.txt",
"FF06955.txt", "FF09419.txt", "FF09476.txt", "FF09533.txt",
"FF09562.txt", "FF09614.txt", "FF09648.txt", "FF09655.txt",
"FF09696.txt", "FF09751.txt", "FF09781.txt", "FF09841.txt",
"FF09901.txt","FF093619.txt","FF096275.txt","FF096413.txt",
"FF096619.txt","FF096688.txt","FF096825.txt","FF096825.txt"
};
vector<string> file_output
{"FF096_out.txt","FF0959_out.txt","FF0971_out.txt","FF0972_out.txt",
"FF06955_out.txt", "FF09419_out.txt", "FF09476_out.txt", "FF09533_out.txt",
"FF09562_out.txt", "FF09614_out.txt", "FF09648_out.txt", "FF09655_out.txt",
"FF09696_out.txt", "FF09751_out.txt", "FF09781_out.txt", "FF09841_out.txt",
"FF09901_out.txt","FF093619_out.txt","FF096275_out.txt","FF096413_out.txt",
"FF096619_out.txt","FF096688_out.txt","FF096825_out.txt","FF096825_out.txt"
};
// tempo, ampiezza, errore
vector <vector <double> > semiampiezze;

if(ForzForzante.size() != file_input.size()) {
    cerr << "Errore: ForzForzante e file_input hanno dimensioni diverse!\n";
    return -1;
}

//string file_input;
/*
 /\_/\
( o.o )
 > ^ <

FOR CHE ANALIZZA TUTTI I DATI
*/

for( int i=0; i< file_input.size(); i++){
    ifstream fin(file_input[i]);
    if (!fin) {
        cerr << "Errore nell'apertura del file di input: " << file_input[i] << endl;
        continue;  // Passa al file successivo invece di terminare
    }

    ofstream fout(file_output[i]);
    if (!fout) {
        cerr << "Errore nell'apertura del file di output: " << file_output[i] << endl;
        continue;
    }




risonanza file;
double temp, forz, pend, amp, fas;
double next_forz = 0;  // Variabile per memorizzare il valore successivo
// Leggi il primo valore fuori dal ciclo
if (fin >> temp >> forz >> pend >> amp >> fas) {
    while (true) {
        // Prova a leggere il valore successivo
        double next_temp, next_pend, next_amp, next_fas;
        if (!(fin >> next_temp >> next_forz >> next_pend >> next_amp >> next_fas)) {
            break;  // Fine del file
        }
          // Controlla sia il valore corrente che il successivo
        if (forz != 0 && next_forz != 0) {
            file.tempo.push_back(temp);
            file.forzante.push_back(forz);
            file.pendolo.push_back(pend);
            file.ampiezza.push_back(amp);
            file.fase.push_back(fas);
        }
        temp = next_temp;
        forz = next_forz;
        pend = next_pend;
        amp = next_amp;
        fas = next_fas;
    }
}


double max_ass = max(file.pendolo);
double min_ass = min(file.pendolo);


// sopra il 80%
vector <int> indici_massimi;
for (int i = 0;  i < file.pendolo.size() ; i++ ){
    if ( file.pendolo[i] >= 0.80*max_ass){
        indici_massimi.push_back(i);
    }
}
// sotto il 80%
vector <int> indici_minimi;
for (int i = 0;  i < file.pendolo.size() ; i++ ){
    if ( file.pendolo[i] <= 0.80*min_ass){
        indici_minimi.push_back(i);
    }
}


vector<pair<double, double>> picchi; // tempo, valore
vector<int> gruppo_corrente;
vector<pair<double, double>> picchi_min; // tempo, valore
vector<int> gruppo_corrente_min;

//massimi
for (size_t i = 0; i < indici_massimi.size(); ++i) {
    gruppo_corrente.push_back(indici_massimi[i]);

    // Se ultimo punto o tempo troppo distante dal successivo
    bool ultimo = (i == indici_massimi.size() - 1);
    bool distanza_grande = false;

    if (!ultimo) {
        double t1 = file.tempo[indici_massimi[i]];
        double t2 = file.tempo[indici_massimi[i + 1]];
        distanza_grande = (t2 - t1 > 0.2);
    }

    if (ultimo || distanza_grande) {
        // Trova massimo nel gruppo
        double max_val = -1e9;
        int max_index = -1;
        for (int idx : gruppo_corrente) {
            if (file.pendolo[idx] > max_val) {
                max_val = file.pendolo[idx];
                max_index = idx;
            }
        }

        if (max_index != -1) {
            double t_max = file.tempo[max_index];
            picchi.push_back({t_max, max_val});
        }

        gruppo_corrente.clear();
    }
}
picchi = dimezza(picchi);
//minimi 
for (size_t i = 0; i < indici_minimi.size(); ++i) {
    gruppo_corrente_min.push_back(indici_minimi[i]);

    // Se ultimo punto o tempo troppo distante dal successivo
    bool ultimo = (i == indici_minimi.size() - 1);
    bool distanza_grande = false;

    if (!ultimo) {
        double t1 = file.tempo[indici_minimi[i]];
        double t2 = file.tempo[indici_minimi[i + 1]];
        distanza_grande = (t2 - t1 > 0.2);
    }

    if (ultimo || distanza_grande) {
        // Trova massimo nel gruppo
        double min_val = 1e9;
        int min_index = -1;
        for (int idx : gruppo_corrente_min) {
            if (file.pendolo[idx] < min_val) {
                min_val = file.pendolo[idx];
                min_index = idx;
            }
        }

        if (min_index != -1) {
            double t_min = file.tempo[min_index];
            picchi_min.push_back({t_min, min_val});
        }

        gruppo_corrente_min.clear();
    }
}
picchi_min = dimezza(picchi_min);

fout << "Dati letti: " << file.tempo.size() << " punti" << endl;
if (file.tempo.empty()) {
    cout << "Nessun dato valido letto!" << endl;
}
fout << "Massimi sopra l'80%: " << indici_massimi.size() << endl;
if (indici_massimi.empty()) {
    cout << "Nessun massimo sopra l'80% trovato!" << endl;
}
fout << "Minimi sotto l'80%: " << indici_massimi.size() << endl;
if (indici_minimi.empty()) {
    cout << "Nessun minimo sotto l'80% trovato!" << endl;
}
fout << "Massimi picchi: " << picchi.size() << endl;
if (picchi.empty()) {
    cout << "Nessun massimo picco!" << endl;
}
fout << "Minimi picchi: " << picchi_min.size() << endl;
if (picchi_min.empty()) {
    cout << "Nessun minimo picco!" << endl;
}fout << "tempo_max\tpendolo_max\ttempo_min\tpendolo_min\n";

// Trova la dimensione massima tra picchi e picchi_min
size_t max_size = max(picchi.size(), picchi_min.size());


double semiampiezza = (mediap(picchi)+abs(mediap(picchi_min)))/2;
double errore_semiampiezza = sqrt(pow((stDevp(picchi)/sqrt(picchi.size())),2)+ pow((stDevp(picchi_min)/sqrt(picchi_min.size())),2));
semiampiezze.push_back({ForzForzante[i], semiampiezza, errore_semiampiezza});
fout << "La media dei picchi massimi e' : " << mediap(picchi)<< endl;
fout << "La deviazione standard dei picchi massimi e' : " << stDevp(picchi)/sqrt(picchi.size()) << endl;
fout << "La media dei picchi minimi e' : " << mediap(picchi_min)<< endl;
fout << "La deviazione standard dei picchi minimi e' : " << stDevp(picchi_min)/sqrt(picchi_min.size()) << endl;
fout << "L'ampiezza delle oscillazioni per la frequenza: " << ForzForzante[i] << " e' " << semiampiezza <<endl;
for (size_t i = 0; i < max_size; ++i) {
    if (i < picchi.size()) {
        fout << picchi[i].first << "\t" << picchi[i].second << "\t";
    } else {
        fout << "N/A\tN/A\t";  // Se non ci sono più massimi
    }

    if (i < picchi_min.size()) {
        fout << picchi_min[i].first << "\t" << picchi_min[i].second;
    } else {
        fout << "N/A\tN/A";  // Se non ci sono più minimi
    }

    fout << "\n";
}

}
/*
 /\_/\
( o.o )
 > ^ <
FOR CHE ANALIZZA TUTTI I DATI
*/
ofstream fout_sem("semiampiezze.txt");
    if (!fout_sem){
        cout << "error output " << endl;
        return -1;
}
for (const auto& row : semiampiezze) {
    for (size_t i = 0; i < row.size(); ++i) {
        fout_sem << row[i];
        if (i != row.size() - 1) fout_sem << "\t";  // Separatore tra colonne
    }
    fout_sem << endl;  // Nuova riga
}
return 0;
}

double media(vector<double> v){
    double S = 0;
    for (auto c : v)
        S += c;
    return S / v.size();
}
double stDev(vector<double> v) {
    double m = media(v);
    double s = 0;
    for (auto c : v) {
        s += (c - m) * (c - m);
    }
    return sqrt(s / (v.size() - 1));
}
double max(vector<double> v){
    double max_val = v[0];
    for (double val : v) {
        if (val > max_val)
            max_val = val;
    }
    return max_val;
}
double min(vector<double> v){
    double min_val = v[0];
        for (double val : v){
            if (val < min_val)
            min_val = val;
    }
return min_val;
}
double somma(const vector<double> v) {
    double sum = 0;
    for (auto c : v) {
        sum += c;
    }
    return sum;
}
vector<pair<double, double>> dimezza(const vector<pair<double, double>>& v) {

    vector<pair<double, double>> dimezzati;
    for (size_t i = 0; i < v.size(); i++) {
        if (i % 2 == 0) {  // Prendi solo gli elementi con indice pari
            dimezzati.push_back(v[i]);
        }
    }
    return dimezzati;
}
double mediap (const vector<pair<double, double>>& v){
 double sum = 0.0;
    for (const auto& p : v) {
        sum += p.second; 
    }
    return sum / v.size();
}
double stDevp(const vector<pair<double, double>>& v) {
    double mean = mediap(v); 
    double som_dif = 0.0;
    
    for (const auto& p : v) {
        som_dif += pow(p.second - mean, 2);
    }
    
    return sqrt(som_dif / (v.size() - 1));
}
