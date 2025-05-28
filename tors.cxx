#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

int main(){

string file_input;
string file_output;

cout << "Nome file input" << endl;
cin >> file_input;

cout << "Nome file output" << endl;
cin >> file_output;

vector <double> tempo ;
vector <double> forzante ;
vector <double> pendolo ;
vector <double> ampiezza ;
vector <double> fase ;
double temp, forz, pend, amp, fas ;

while (fin >> temp >> forz >> pend >> amp >> fas){
	tempo.push_back(temp);
	forzante.push_back(forz);
	pendolo.push_back(pend);
	ampiezza.push_back(amp);
	fase.push_back(fas);
}

double max =0;
for (auto c : pendolo){
if (c > max){
max = c;
}}
// sopra il 80%



return 0;
}
