#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <cmath>
using namespace std;

struct Data{
    int seq_num;
    int pos;
    string olignuc;

};

struct Graph{
    int ver;
    vector<vector<int> > adj_list;
    vector<int> n_amount;

    Graph(int vertices): ver(vertices), adj_list(vertices), n_amount(vertices){};

    void addVertex(int first, int second)
    {
    this->adj_list[first].push_back(second);
    this->adj_list[second].push_back(first);
    n_amount[first]++;
    n_amount[second]++;
    }
};

vector<string> sequences;
vector<string> old_sequences;
vector<vector<int> > quali;
vector<string> id;
vector<int> sq_num(5);
vector<vector<pair<char, int> > >old_pos(5);
vector<Data> v_info;
set<int> maxCliqueVertices;
int maxx =0;

void INPUT();
void CliqueHeu(Graph& G, vector<int>& U, int sizee, set<int>& currentClique );
void Remove(int level);
bool LengthCheck();
void Graph_maker(int len,Graph&  G);
void MaxCliqeHeu(Graph& G);

int main(){

     INPUT();

    for(int i=0;i<5;i++){
        cout<<"id: "<<i<<" "<<id[i]<<endl;
        cout<<"sequence: "<<i<<" "<<sequences[i]<<endl;
        cout << "quality: " << i << " ";
        for (int j = 0; j < quali[i].size(); j++) {
            cout << quali[i][j] << " ";
        }
        cout<<endl;
    }

    if(LengthCheck()){
        cout<<"Data entered successfully"<<endl;
     }
     else{
        cout<<"Error when downloading data"<<endl;
        return 0;
     }

     int level;
     cout<<"Enter the quality value: "<<endl;
     cin>>level;
     Remove(level);
     for(int i=0; i<5; i++){
        cout<<"seq after remove: "<<sequences[i]<<endl;
     }
     for(int i=0; i<5; i++){
        for(int j=0; j<sequences[i].size(); j++){
            cout<<"("<<old_pos[i][j].first<<","<<old_pos[i][j].second<<")";
        }
        cout<<endl;
    }

    int len, ver =0;
    cout<<"Enter the length of the oligonucleotides: "<<endl;
    cin>>len;
    for(int i=0; i<5; i++){
       ver+=sequences[i].size()-len+1;
    }
    cout<<"amount vertex: "<<ver<<endl;
    Graph G(ver);
    Graph_maker(len, G);


    cout << "Graph vertices: " << G.ver << endl;
    for (int i = 0; i < G.ver; ++i) {
        cout << "Adjacent vertices to vertex " << i << ": ";
        for (int j=0; j<G.n_amount[i]; j++) {
            cout << G.adj_list[i][j] << " ";
        }
        cout << endl;
    }

    cout << "v_info:" << endl;
    for (size_t i = 0; i < v_info.size(); ++i) {
        cout <<i<<": "<<"seq_num: " << v_info[i].seq_num << ", pos: " << v_info[i].pos << ", olignuc: " << v_info[i].olignuc << endl;
    }

    MaxCliqeHeu(G);
    cout<<"max clique size: "<<maxx<<endl;
    cout<<"max clique: "<<endl;
     for (int vertex : maxCliqueVertices) {
        cout << vertex << ": "<<" seq: "<<v_info[vertex].seq_num<<" pos: "<<v_info[vertex].pos<<" olic: "<<v_info[vertex].olignuc<<endl;
    }
}

void INPUT() {
     string line;
    int line_num = 0;

    ifstream plik("fasta.txt");
    if (plik.is_open()) {
            line_num = 0;

    while (getline(plik, line)) {
                line_num++;

            //cout << "wykryto linie nr: " << line_num<<"linia:"<< line << endl;
            size_t space_pos = line.find(' ');
            string fragment = line.substr(0,space_pos);
            //cout << fragment << endl;
            id.push_back(fragment);
            size_t start_pos = line.find("38_ ");
            string fragment1 = line.substr(start_pos + 4);
            sequences.push_back(fragment1);
}


            plik.clear();
            plik.seekg(0, ios::beg);
        }
        plik.close();

    ifstream plik1("qual.txt");
    if (plik1.is_open()) {
    for(int i =0; i<5; i++){
        line_num = 0;
        string line;
        getline(plik1, line);
        istringstream iss(line);
        int value;
        vector<int> row;
        while(iss>>value){
            row.push_back(value);
        }
        quali.push_back(row);
    }
    }

    plik1.close();
    }

void MaxCliqeHeu(Graph& G){
    for(int i=0; i< G.ver; i++){
        if(G.n_amount[i] >= maxx){
            vector<int> U;
            for(int j=0; j<G.adj_list[i].size(); j++){
                if(G.n_amount[G.adj_list[i][j]]>=maxx){
                    U.push_back(G.adj_list[i][j]);
                }
            }
            set<int> currentClique;
            currentClique.insert(i);
            CliqueHeu(G,U,1, currentClique);
        }
    }
}

void CliqueHeu(Graph& G, vector<int>& U, int sizee, set<int>& currentClique){
    if(U.empty()){
        if(sizee>maxx){
            maxx=sizee;
            maxCliqueVertices = currentClique;
        }
        return;
    }

    int maxdegree = -1;
    int maxdegree_index = -1;

    for(int i=0; i<U.size(); i++){
       if(G.n_amount[U[i]]>maxdegree) {
        maxdegree=G.n_amount[U[i]];
        maxdegree_index = i;
       }
    }
    int u = U[maxdegree_index];
    U.erase(remove(U.begin(), U.end(), u), U.end());
    vector<int> N_prime;
    for (int w : G.adj_list[u]) {
        if (G.n_amount[w] >= maxx) {
            N_prime.push_back(w);
        }
    }
    currentClique.insert(u);
    CliqueHeu(G, U, sizee + 1, currentClique);
    currentClique.erase(u);
}

void Remove(int level){
    for(int i=0; i<5; i++){
        for(int j=0; j<sequences[i].size(); j++){
            old_pos[i].push_back(make_pair(sequences[i][j], j+1));
        }
    }
    for(int i=0; i<5; i++){
        for(int j=0; j<sequences[i].size(); j++){
            cout<<"("<<old_pos[i][j].first<<","<<old_pos[i][j].second<<")";
        }
        cout<<endl;
    }

    old_sequences = sequences;

    for(int i=0; i<5; i++){
            int a=0;
        for(int j=0; j<old_sequences[i].size(); j++){
            if(quali[i][j]<level){
                old_pos[i].erase(old_pos[i].begin() + a);
                sequences[i].erase(a,1);
                a--;
            }
            a++;
        }
    }
    cout<<"done"<<endl;
}

bool LengthCheck(){
    for(int i=0; i<5; i++){
        if(sequences[i].size()!=quali[i].size()){
            cout<<i<<": "<<sequences[i].size()<<", "<<quali[i].size()<<endl;
            return 0;
        }
    }
    return 1;
}

void Graph_maker(int len,Graph&  G){
    cout << "Inside Graph_maker function" << endl;
    for(int i=0; i<5; i++){
        for(int j=0; j<sequences[i].size()-len+1; j++){
            string seq;
            int seqn;
            int pos;
            seq = sequences[i].substr(j,len);
            seqn = i;
            pos = old_pos[i][j].second;
            v_info.push_back({seqn, pos, seq});
        }
    }
    for (int i = 0; i < v_info.size(); i++) {
        for (int j = 0; j < v_info.size(); j++) {
            if(i!=j){
            int pos1 = v_info[i].pos;
            int pos2 = v_info[j].pos;
            int diff = abs(pos1 - pos2);
        if (v_info[i].olignuc == v_info[j].olignuc && v_info[i].seq_num != v_info[j].seq_num && diff<=len*10){
            cout << "Adding vertex (" << i << ", " << j << ") to the graph" << endl;
            G.addVertex(i, j);
            }
        }
    }
}
    cout << "Exiting Graph_maker function" << endl;

}

