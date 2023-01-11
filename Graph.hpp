#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <vector>
#include <list>
#include <queue>
#include <limits.h>
#include<string>
#include <iostream>
#include <bits/stdc++.h>

using namespace std;

template <class T> class Edge;
template <class T> class Graph;

  
struct lista{
	int from;
	int to;
	lista(int from, int to)
    {
        this->from= from;
		this->to= to;
    }
};

typedef  pair<int, int> iPair;
  

struct praph{
	
    int V, E;
    vector< pair<int, iPair> > edges;
	vector<lista> akmes;
    praph(int V, int E)
    {
        this->V = V;
        this->E = E;
    }
  
	void addAkmi(int from, int to){
		akmes.push_back(lista(from, to));
	}
    
    void addEdge(int u, int v, int w)
    {
        edges.push_back({w, {u, v}});
    }
  
    
    void kruskalMST();
};
  

struct DisjointSets
{
    int *parent, *rnk;
    int n;
  
    
    DisjointSets(int n)
    {
        
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];
  
        
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;
  
            
            parent[i] = i;
        }
    }
  
    
    int find(int u)
    {
        
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }
  
    
    void merge(int x, int y)
    {
        x = find(x), y = find(y);
  
        
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else 
            parent[x] = y;
  
        if (rnk[x] == rnk[y])
            rnk[y]++;
    }
};

template <typename T>
class Vertex {
    bool addEdge(Vertex<T> *from, Vertex<T> *to, int w);
	bool removeEdgeTo(Vertex<T> *d);
	
public:
	/*~Vertex(){
		typename vector<Edge<T> >::iterator it= adj.begin();
		typename vector<Edge<T> >::iterator ite= adj.end();
		for (; it!=ite; it++){
			//delete &it;
		}
	}*/
	vector<Edge<T>> adj;
	Vertex(T in);
    bool visited;
	int row;
    T info;
    T getInfo() const;
    friend class Graph<T>;
};

template <class T>
Vertex<T>::Vertex(T in): info(in){}

template <class T>
T Vertex<T>::getInfo() const {
	return info;
}

template<typename T>
struct Edge {
  Vertex<T> *from;
  Vertex<T> *to;
  int dist;
  int row;
  Edge(Vertex<T> *f, Vertex<T> *t, int d): from(f), to(t), dist(d){};
  bool operator<(const Edge<T>& e) const;
  bool operator>(const Edge<T>& e) const;

  template<typename U>
  friend std::ostream& operator<<(std::ostream& out, const Edge<U>& e);
};

template<typename T>
std::ostream& operator<<(std::ostream& out, const Edge<T>& e) {
  out << e.from << " -- " << e.to << " (" << e.dist << ")";
  return out;
}

template <class T>
bool Vertex<T>::removeEdgeTo(Vertex<T> *d) {
	typename vector<Edge<T> >::iterator it= adj.begin();
	typename vector<Edge<T> >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->from->info == d->info) {
			adj.erase(it);
			return true;
		}
		else it++;
	}
	return false;
}


template <class T>
bool Vertex<T>::addEdge(Vertex<T> *from, Vertex<T> *to, int w) {

    typename vector<Edge<T>>::iterator it= adj.begin();
	typename vector<Edge<T>>::iterator ite= adj.end();
	
	int found=0;
	while (found!=1 && it!=ite ) {
		if ( it->from->info == from->info && it->to->info== to->info ){
            found++;
        }
		it ++;
	}
	if(found!=0){
        return(false);
    }
    
	Edge<T> edgeD(from, to, w);
	edgeD.row= edgeD.to->row;
	adj.push_back(edgeD);
	
	std::sort(adj.begin(), adj.end(),[](const Edge<T>& lhs, const Edge<T>& rhs){
		return lhs.row < rhs.row;
	});
    return(true);
}


template <typename T>

class Graph {

vector<Vertex<T> *> arrayVertex;
bool directed;
public:
  Graph(bool isDirectedGraph){
	  if(isDirectedGraph== true){
		  directed= true;
	  }
	  else{
		  directed= false;
	  }
  }
  ~Graph(){
	typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::iterator ite= arrayVertex.end();
    for (; it!=ite; it++){
		delete *it;
	}
	//delete *arrayVertex;
  }
  bool contains(const T& info);
  bool addVtx(const T& info);
  bool rmvVtx(const T& info);
  bool addEdg(const T& from, const T& to, int distance);
  bool rmvEdg(const T& from, const T& to);
  std::list<Vertex<T>> dfs(const T& info);
  std::list<Vertex<T>> bfs(const T& info) const;
  std::list<Edge<T>> mst();

  bool print2DotFile(const char *filename) const;
  std::list<Vertex<T>> dijkstra(const T& from, const T& to);
};


template <typename T>
bool Graph<T>::contains(const T& info){

    typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::iterator ite= arrayVertex.end();
    for (; it!=ite; it++){
        if((*it)->info== info){
            return(true);
        }
    }
    return(false);
}

template <typename T>
bool Graph<T>::addVtx(const T& info){

    if(contains(info)== false){
        Vertex<T> *v1 = new Vertex<T>(info);
		v1->row= arrayVertex.size();
        arrayVertex.push_back(v1);
        return(true);
    }
    else{
        return(false);
    }
}

template <typename T>
bool Graph<T>::rmvVtx(const T& info){

    if(contains(info)== true){

        typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
        while((*it)->info!= info){
            it++;
        }
        Vertex<T> * v= *it;
        arrayVertex.erase(it);
        typename vector<Vertex<T>*>::iterator it1= arrayVertex.begin();
        typename vector<Vertex<T>*>::iterator it1e= arrayVertex.end();
        for (; it1!=it1e; it1++) {
            (*it1)->removeEdgeTo(v);
        }
        delete v;
        return true;
    }
    return(false);
}

template <class T>
bool Graph<T>::addEdg(const T& from, const T& to, int cost){

    typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::iterator ite= arrayVertex.end();
	int found=0;
	Vertex<T> *edgeS, *edgeD;
	
	while (found!=2 && it!=ite ) {
		
		if ( (*it)->info == from )
			{ edgeS=*it; found++;}
		if ( (*it)->info == to )
			{ edgeD=*it; found++;}
		it ++;
	}
	if (found !=2){
		return false;
	}
	if(directed== false){
		bool hey= edgeD->addEdge(edgeD,edgeS,cost);
		if(hey== true){
		}
	}
	return edgeS->addEdge(edgeS,edgeD,cost);
}

template <typename T>
bool Graph<T>::rmvEdg(const T& from, const T& to){

    typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::iterator ite= arrayVertex.end();
	int found=0;
	Vertex<T> *edgeS, *edgeD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == from )
			{ edgeS=*it; found++;}
		if ( (*it)->info == to )
			{ edgeD=*it; found++;}
		it ++;
	}
	if (found!=2) return false;
	
	if(directed== false){
		bool hey=edgeD->removeEdgeTo(edgeS);
		if(hey==false){
		}
	}
	return edgeS->removeEdgeTo(edgeD);
}

template <typename T>
std::list<Vertex<T>> Graph<T>::dfs(const T& info) {
	typename vector<Vertex<T>*>::const_iterator it= arrayVertex.begin();// or const iterator
	typename vector<Vertex<T>*>::const_iterator ite= arrayVertex.end();
	for (; it !=ite; it++){
		(*it)->visited=false;
    }
    
    it= arrayVertex.begin();
    int j=9;
    for (; j!=10; it++){
        if((*it)->info== info){
            j=10;
        }
    }
	it--;
	list<Vertex<T>> res;
	for (; it !=ite; it++){
		if ( (*it)->visited==false ){
			dfs2(*it,res);
		}
	}
	return res;
}

template <typename T>
void dfs2(Vertex<T> *v, list<Vertex<T>> &res) {// ή ....const{
	
    v->visited = true;
	res.push_back(*v);
	typename vector<Edge<T> >::iterator it= (v->adj).begin();
	typename vector<Edge<T> >::iterator ite= (v->adj).end();
	for (; it !=ite; it++){
	    if ( it->to->visited == false ){
			dfs2(it->to, res);
        }
    }
}

template <typename T>
std::list<Vertex<T>> Graph<T>::bfs(const T& info) const{
	
	typename vector<Vertex<T>*>::const_iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::const_iterator ite= arrayVertex.end();
	Vertex<T> *k;
	Vertex<T> *counter;
	while ( it!=ite ) {
		counter= *it;
		if ( (*it)->info == info ){
			k=*it;
		}
		counter->visited= false;
		it ++;
	}
	
	list<Vertex<T>> res;
	queue<Vertex<T> *> q;
	q.push(k);
	k->visited = true;
	
	while (!q.empty()) {
		Vertex<T> *v1 = q.front();
		q.pop();
		res.push_back(v1->info);
		typename vector<Edge<T> >::iterator it2=v1->adj.begin();
		typename vector<Edge<T> >::iterator ite2=v1->adj.end();
		
		for (; it2!=ite2; it2++) {
			Vertex<T> *d = it2->to;
			if (d->visited==false) {
				d->visited=true;
				q.push(d);
			}
		}
	}
	return res;
}

template <typename T>
std::list<Edge<T>> Graph<T>::mst(){

    list<Edge<T>> res;
    vector<Edge<T>> edges;
	int j=5;
	if(directed== true){
		return(res);
	}
    typename vector<Vertex<T>*>::const_iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::const_iterator ite= arrayVertex.end();
	
	
	
    for (; it !=ite; it++){
		
        typename vector<Edge<T> >::iterator it1= ((*it)->adj).begin();
        typename vector<Edge<T> >::iterator ite1= ((*it)->adj).end();
		
        for (; it1 !=ite1; it1++){
			typename vector<Edge<T> >::iterator it2= edges.begin();
			typename vector<Edge<T> >::iterator ite2= edges.end();
			
			for (; it2 !=ite2; it2++){
				if(it2->from->info== it1->to->info && it2->to->info== it1->from->info){
					j=10;
				}
			}
			if(j!=10){
				Edge<T> edgeD(it1->from, it1->to, it1->dist);
				edgeD.row= edges.size();
				edges.push_back(edgeD);
			}
			j=5;
        }
    }
    int V = arrayVertex.size(), E = edges.size();
	praph g(V, E);
  
	typename vector<Edge<T> >::iterator it2= edges.begin();
	typename vector<Edge<T> >::iterator ite2= edges.end();
    for(; it2!= ite2; it2++){
		g.addEdge(it2->from->row,it2->to->row, it2->dist);
	}
  
    g.kruskalMST();
	
	typename vector<lista>::iterator ip= g.akmes.begin();
	typename vector<lista>::iterator ipe= g.akmes.end();
	for(; ip!=ipe; ip++){
		it2= edges.begin();
		for(; it2!= ite2; it2++){
			Edge<T> akmi(it2->from, it2->to, it2->dist);
			//typename vector<lista >::iterator test= ip++;
			
			int from= ip->from;
			int to= ip->to;
			if(from== it2->from->row && to== it2->to->row){
				akmi.from= it2->from;
				akmi.to= it2->to;
				akmi.dist= it2->dist;
				res.push_back(akmi);
			}
		}
	}
		
    return(res);
}

void praph::kruskalMST()
{
    int mst_wt = 0; 
    sort(edges.begin(), edges.end());
	vector<lista> g1;
    
    DisjointSets ds(V);
  
    vector< pair<int, iPair> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++){
		
        int u = it->second.first;
        int v = it->second.second;
  
        int set_u = ds.find(u);
        int set_v = ds.find(v);
  
        
        if (set_u != set_v)
        {
			lista edgeD(u, v);
            akmes.push_back(edgeD);
  
            
            mst_wt += it->first;
  
            ds.merge(set_u, set_v);
        }
    }
}
  


template <typename T>
int minDistance(int dist[],bool sptSet[], vector<Vertex<T> *> &arrayVertex){

    int min = INT_MAX;
    int min_index;
    int V= arrayVertex.size();

    for (int v = 0; v < V; v++){
        if (sptSet[v] == false && dist[v] <= min){
            min = dist[v];
            min_index = v;
        }
    }
    return min_index;
}

template <typename T>
void Path(int parent[], int j, list<Vertex<T>> &res, vector<Vertex<T> *> &arrayVertex)
{

    if (parent[j] == - 1)
        return;

    Path(parent, parent[j], res, arrayVertex);

    typename vector<Vertex<T>*>::iterator ite= arrayVertex.begin();
    for(int i=0; i< j && ite!= arrayVertex.end(); i++){
        ite++;
    }
    Vertex<T> *v= *ite;
	
    res.push_back(*v);
}

template <typename T>
list<Vertex<T>> Graph<T>::dijkstra( const T& from, const T& to){

    int V= arrayVertex.size();
    int graph2[V][V];
    int dist[V];
    list<Vertex<T>> res;
    bool sptSet[V];
    int parent[V];

    for(int i=0; i<V; i++){
        for(int j=0; j<V; j++){
            graph2[i][j]=0;
        }
    }

    typename vector<Vertex<T>*>::iterator it= arrayVertex.begin();
	typename vector<Vertex<T>*>::iterator ite= arrayVertex.begin();

    for(int i=0; i<V && it!= arrayVertex.end(); i++){

        typename vector<Edge<T>>::iterator it1= ((*it)->adj).begin();
        typename vector<Edge<T>>::iterator ite1= ((*it)->adj).end();

        for (; it1 !=ite1; it1++){
            int f=0;
            Vertex<T> *k= it1->to;
            while((*ite)->info!= k->info){
                f++;
                ite++;
            }
            ite= arrayVertex.begin();

            graph2[i][f]= it1->dist;
            f=0;
        }
        it++;
    }
    
   /* for(int i=0; i<V; i++){
		for(int j=0; j<V; j++){
			cout<< graph2[i][j]<< " ";
		}
		cout<< "\n";
	}*/


    for (int i = 0; i < V; i++){
        parent[i] = -1;
        dist[i] = INT_MAX;
        sptSet[i] = false;
    }

    ite= arrayVertex.begin();
    int f=0;
    while((*ite)->info!= from){
        f++;
        ite++;
    }
    dist[f] = 0;

    for (int count = 0; count < V - 1; count++){

        int u = minDistance(dist, sptSet, arrayVertex);
        sptSet[u] = true;
        for (int v = 0; v < V; v++){
            if (!sptSet[v] && graph2[u][v] && dist[u] + graph2[u][v] < dist[v]){
                parent[v] = u;
                dist[v] = dist[u] + graph2[u][v];
            }
        }
    }

    ite= arrayVertex.begin();
    int j=0;
    while((*ite)->info!= to){
        j++;
        ite++;
    }

    if (parent[j] == - 1 || dist[j]== INT_MAX){
        return res;
    }

    ite= arrayVertex.begin();
    while((*ite)->info!= from){
        ite++;
    }
	Vertex<T> *v= *ite;
    res.push_back(*v);
	
    Path(parent, j, res, arrayVertex);
	/*for(int i=0; i<V; i++){
		cout<< parent[i];
		cout<<" ";
	}*/
    return(res);
}


template <typename T>
bool Graph<T>::print2DotFile(const char *filename) const{

    string bigstring= ("graph Tree {\n ");
    int V= arrayVertex.size();

    typename vector<Vertex<T>*>::const_iterator it= arrayVertex.begin();

    for(int i=0; i<V; i++){
        typename vector<Edge<T>>::const_iterator it1= ((*it)->adj).begin();
        typename vector<Edge<T>>::const_iterator ite1= ((*it)->adj).end();
        for(; it1!=ite1; it1++){
            std::cout << bigstring << "\t";
            Vertex<T> *j= it1->from;
            std::cout << j->getInfo();
            std::cout << " -- ";
            Vertex<T> *k= it1->to;
            std::cout << k->getInfo();
            std::cout << "\n";
        }
        it++;
    }

    std::cout << "}";
    std::cout << "\n";
    return (true);
}

#endif
/*
clear
g++ -Wall -g -std=c++11 GraphString.cpp -o gstr
g++ -Wall -g -std=c++11 GraphStudent.cpp -o gstudent
g++ -Wall -g -std=c++11 GraphInteger.cpp -o gint
./gstr < tests/dijkstra1.in > dijkstra1.out
diff -urN tests/dijkstra1.std dijkstra1.out
*/
