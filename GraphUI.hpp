
#include <iostream>
#include <string>
#include <cstring>
#ifndef _GRAPH_UI_
#define _GRAPH_UI_

template <typename T>
int graphUI() {

  string option, line;
  int distance;
  bool digraph = false;

  cin >> option;
  if(!option.compare("digraph"))
    digraph = true;
  Graph<T> g(digraph);

  while(true) {

    std::stringstream stream;
    cin >> option;

    if(!option.compare("av")) {
      getline(std::cin, line);
      stream << line;
      T vtx(stream);
      if(g.addVtx(vtx))
        cout << "av " << vtx << " OK\n";
      else
        cout << "av " << vtx << " NOK\n";
    }
    else if(!option.compare("rv")) {
        getline(std::cin, line);
        stream << line;
        T vtx(stream);
        if(g.rmvVtx(vtx))
            cout << "rv " << vtx << " OK\n";
        else
            cout << "rv " << vtx << " NOK\n";
    }
    else if(!option.compare("ae")) {
        getline(std::cin, line);
        stream << line;
        T vtx(stream);
        T vtx1(stream);
        T dis(stream);

        stringstream geek;
        geek << dis;

        distance = 0;
        geek >> distance;
		
        if(g.addEdg(vtx, vtx1, distance))
            cout << "ae " << vtx << " " << vtx1 << " OK\n";
        else
            cout << "ae " << vtx << " " << vtx1 << " NOK\n";

    }
    else if(!option.compare("re")) {
        getline(std::cin, line);
        stream << line;
        T vtx(stream);
        T vtx1(stream);

        if(g.rmvEdg(vtx, vtx1))
            cout << "re " << vtx << " " << vtx1 << " OK\n";
        else
            cout << "re " << vtx << " " << vtx1 << " NOK\n";

    }
    else if(!option.compare("dot")) {

        if(g.print2DotFile("C:\\Users\\dimit\\hw5_skeleton\\dot"))
            cout << "dot file OK\n";
        else
            cout << "dot file NOK\n";

    }
    else if(!option.compare("bfs")) {
        getline(std::cin, line);
        stream << line;
        T vtx(stream);

        list <Vertex<T>> res= g.bfs(vtx);

        cout << "\n----- BFS Traversal -----\n";

        typename list<Vertex<T>>::iterator it= res.begin();
        typename list<Vertex<T>>::iterator ite= res.end();
        cout << it->info;
        if(res.size()>1){
            it++;
            for (; it !=ite; it++){
                cout << " -> " << it->info;
            }
        }

        cout << "\n-------------------------\n";
    }
    else if(!option.compare("dfs")) {
        getline(std::cin, line);
        stream << line;
        T vtx(stream);

        list<Vertex<T>> res= g.dfs(vtx);

        cout << "\n----- DFS Traversal -----\n";
        typename list<Vertex<T>>::iterator it= res.begin();
        typename list<Vertex<T>>::iterator ite= res.end();
        cout << it->info;
        if(res.size()>1){
            it++;
            for (; it !=ite; it++){
                cout << " -> " << it->info;
            }
        }

        cout << "\n-------------------------\n";
    }
    else if(!option.compare("dijkstra")) {
      getline(std::cin, line);
      stream << line;
      T from(stream);
      T to(stream);

      list<Vertex<T>> res= g.dijkstra(from, to);
      cout << "Dijkstra (" << from << " - " << to <<"): ";

      typename list<Vertex<T>>::iterator it= res.begin();
      typename list<Vertex<T>>::iterator ite= res.end();
      cout << it->info;
      if(res.size()>1){
        it++;
        for (; it !=ite; it++){
            cout << ", " << it->info;
        }
      }
      cout<<"\n";

    }
    else if(!option.compare("mst")) {

        list<Edge<T>> res= g.mst();
        cout << "\n--- Min Spanning Tree ---\n";
        typename list<Edge<T>>::iterator it= res.begin();
        typename list<Edge<T>>::iterator ite= res.end();
        int sum=0;
        for (; it !=ite; it++){
			sum= sum+ it->dist;
            cout << it->from->info << " -- " << it->to->info << " (" << it->dist << ")";
			cout << "\n";
        }

        cout << "MST Cost: " << sum << endl;
    }
    else if(!option.compare("q")) {
      cerr << "bye bye...\n";
      return 0;
    }
    else if(!option.compare("#")) {
      string line;
      getline(cin,line);
      cerr << "Skipping line: " << line << endl;
    }
    else {
      cout << "INPUT ERROR\n";
      return -1;
    }
  }
  return -1;
}

#endif
