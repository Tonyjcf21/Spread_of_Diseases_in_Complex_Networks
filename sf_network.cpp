#include<iostream>
#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<time.h>
#include<ctime>
#include<random>
#include<thread>
#include<sstream>
#include <algorithm>
//----
#include <vector>
#include <list>
#include <set>
#include <iterator>
//----
//#include "gnuplot.h"

using namespace std;

#define TAU 6.283185307

// g++ adj_list44_novisitados.cpp -o sfn2 -std=c++2a

float beta;
int betastr;

class Edge;
class Node;
void tostr(float);

//OJO que X e Y deben ser la longitud del edge, no?
//podrias definir el vector distancia r
//no se si nos interesan las distancias aqui.!!!!!!!
class Edge
{
	public:
		int destinationID;
		float x,y; //position
		int n_edges = 0;

		Edge() {}
		Edge(int destID){
			destinationID = destID;
		}

		//void setEdgeValues(int destID){
		//	destinationID = destID;
		//}

		void setPosition(float posX, float posY){
			x = posX;
			y = posY;
			//sqrt((P[i].x-P[j].x)*(P[i].x-P[j].x) + (P[i].y-P[j].y)*(P[i].y-P[j].y)
		}

		int getDestinationID(){
			return destinationID;
		}

		/*
		float getx(){
			return x;
		}

		float gety(){
			return y;
		}
		*/
};

//CAPAZ PONER EL STATUS AQUI EN Node???
//SI, EL STATUS ES DE LOS NODOS O NodeA, NO DEL EDGE.

class Node {

public:
	int node_id;
	int status; //0=S, 1=I, 2=D
	int recovery_time;
	int dead_time;
	list <Edge> edgeList;

	Node(){
		node_id = 0;
		status = 0;
		recovery_time = 10;
		dead_time = 2;
	}

	/*
	Node(int id, int Stat, int rec_t, int dead_t){
		node_id = id;
		status = Stat;
		recovery_time = rec_t;
		dead_time = dead_t;
	}
	

	void setStatus(){
		status = Stat;
	}
	*/
	void printEdgeList() {
    	cout << "[";
    	for (auto it = edgeList.begin(); it != edgeList.end(); it++) {
      		cout << it -> getDestinationID() << " --> ";
    	}
    	cout << "]";
    	cout << endl;
    }

    void printNumberEdges(){
    	cout << "Number of Edges for this node: " << edgeList.size() << endl;
    }


	void setID(int id){
		node_id = id;
	}

	int getID(){
		return node_id;
	}

	list <Edge> getEdgeList() {
		return edgeList;
	}

	void addEdgeToEdgelist(int toNodeID)
    {
  		Edge e(toNodeID);
  		edgeList.push_back(e); 
  		//cout<<"Edge between "<<state_id<<" and "<<toNodeID<<" added Successfully"<<endl; 	
    }
/* Weighted network. I dont have getWeight() in Edge class.

    void printEdgeList() {
    cout << "[";
    	for (auto it = edgeList.begin(); it != edgeList.end(); it++) {
    		cout << it -> getDestinationID() << "(" << it -> getWeight() << ") --> ";
    	}
    cout << "]";
    cout << endl;
    }
*/
    void clearEdges(){
    	edgeList.clear();
    }
};

class Graph{

	vector<Node> nodes;

public:
	//vector<Node> nodes;

	void addNode(Node newNode){
		bool check = checkIfNodeExistsByID(newNode.getID());
		if (check == true)
		{
			cout << "Node already exists" << endl;
			//break;
		}
		else
		{
			nodes.push_back(newNode);
		}
	}

	int NumberOfEdgesByNode(int currentNode){
		int temp2;

		temp2 = nodes.at(currentNode).edgeList.size();
		return temp2;
	}

	float DegreeOfNodeToAttach(int node_to_attach){
		return nodes.at(node_to_attach).edgeList.size();
	}

	//Removed weight from addEdgeID.
	void addEdgeByID(int fromNode, int toNode) {
		int CurrentEdgeNumber = 0;
    	bool check1 = checkIfNodeExistsByID(fromNode);
    	bool check2 = checkIfNodeExistsByID(toNode);

    	bool check3 = checkIfEdgeExistsByID(fromNode, toNode);
    	if ((check1 && check2 == true)) {

        	if (check3 == true) {
        		cout << "Edge between " << getNodeByID(fromNode).getID() << "(" << fromNode << ") and " << getNodeByID(toNode).getID() << "(" << toNode << ") Already Exist" << endl;
      		} 
      		else {
				for (int i = 0; i < nodes.size(); i++) {

            		if (nodes.at(i).getID() == fromNode) {
            			//cout << "node == from... ojo aqui " << "from: " << fromNode << "... to: " << toNode << endl;
            			Edge e(toNode); //Edge e(toNode, weight);
            			Edge e2(fromNode);
            			//edgeList.push_back(e); 
            			//nodes.at(i).addEdgeToEdgelist(toNode,weight);
            			nodes.at(i).edgeList.push_back(e);
            			CurrentEdgeNumber = nodes.at(i).edgeList.size();
            			nodes.at(toNode).edgeList.push_back(e2);
          			} 
          			/*
          			else if (nodes.at(i).getID() == toNode) {
          				cout << "node == to... ojo aqui " << "from: " << fromNode << "... to: " << toNode << endl;
            			Edge e(toNode); //Edge e(toNode, weight);
		            	//edgeList.push_back(e); 
		            	//nodes.at(i).addEdgeToEdgelist(fromNode,weight);
		            	nodes.at(i).edgeList.push_back(e);
            		}
            		*/
        		}
        		//cout << "Edge between " << fromNode << " and " << toNode << " added Successfully" << endl;
        	}
    	} 
    	else {
    		cout << "Invalid Node ID entered.";
    	}
    }

	bool checkIfNodeExistsByID(int nID){
		bool flag = false;

		for (int i = 0; i < nodes.size(); ++i)
		{
			if (nodes.at(i).getID()==nID)
			{
				return true;
			}
		}
		return flag;
	}

	bool checkIfEdgeExistsByID(int fromNode, int toNode) {
    	Node v = getNodeByID(fromNode);
    	list < Edge > e;
    	e = v.getEdgeList();
    	bool flag = false;
    	
    	for (auto it = e.begin(); it != e.end(); it++) {
      		if (it -> getDestinationID() == toNode) {
        		flag = true;
        		return flag;
        		break;
      		}

    	}

    	Node v2 = getNodeByID(toNode);
    	list < Edge > e2;
    	e2 = v2.getEdgeList();
    	
    	for (auto it = e2.begin(); it != e2.end(); it++) {
      		if (it -> getDestinationID() == fromNode) {
        		flag = true;
        		return flag;
        		break;
      		}

    	}
    	return flag;
  	}


	Node getNodeByID(int vid) {
		Node temp;
		for (int i = 0; i < nodes.size(); i++) {
			temp = nodes.at(i);
			if (temp.getID() == vid) {
			return temp;
			}
		}
		return temp;
	}

	void clearGraph(){
		for (int i = 0; i < nodes.size(); i++) {
		    Node temp;
		    temp = nodes.at(i);
		    //cout << nodes.size();
		    //cout << temp << endl;
		    //cout << " (" << temp.getID() << ") --> ";
		    temp.clearEdges();
		    
		}
		nodes.clear();
	}

	void printGraph() {
	    for (int i = 0; i < nodes.size(); i++) {
		    Node temp;
		    temp = nodes.at(i);
		    //cout << nodes.size();
		    //cout << temp << endl;
		    cout << " (" << temp.getID() << ") --> ";
		    temp.printEdgeList();
		    temp.printNumberEdges();
		    /*
		    cout << "[";
    	for (auto it = edgeList.begin(); it != edgeList.end(); it++) {
      		cout << it -> getDestinationID() << " --> ";
    	}
    	cout << "]";
    	cout << endl;

    void showlist() //(list<int> g)
	{
    	list<int>::iterator it;

    	for (auto it = edgeList.begin(); it != edgeList.end(); ++it)
        	cout << '\t' << *it;

    	cout << '\n';
	}
    	*/
	    }
    }

	//
	/*
	void updateEdgeByID(int fromVertex, int toVertex, int newWeight) {
    bool check = checkIfEdgeExistsByID(fromVertex, toVertex);
    if (check == true) {
      for (int i = 0; i < vertices.size(); i++) {

        if (vertices.at(i).getStateID() == fromVertex) {
          for (auto it = vertices.at(i).edgeList.begin(); it != vertices.at(i).edgeList.end(); it++) {
            if (it -> getDestinationVertexID() == toVertex) {
              it -> setWeight(newWeight);
              break;
            }

          }

        } else if (vertices.at(i).getStateID() == toVertex) {
          for (auto it = vertices.at(i).edgeList.begin(); it != vertices.at(i).edgeList.end(); it++) {
            if (it -> getDestinationVertexID() == fromVertex) {
              it -> setWeight(newWeight);
              break;
            }

          }
        }
      }
      cout << "Edge Weight Updated Successfully " << endl;
    } else {
      cout << "Edge between " << getVertexByID(fromVertex).getStateName() << "(" << fromVertex << ") and " << getVertexByID(toVertex).getStateName() << "(" << toVertex << ") DOES NOT Exist" << endl;
    }
  }
*/

};

/*
The Algorithm:
Input: Number of Nodes N; 
       Initial number of nodes m0; 
       Offset Exponent a; 
       Minimum degree 1 <= d <= m0.
Output: scale-free multigraph G = ({0,....,N-1}, E).

1) Add m0 nodes to G.
2) Connect every node in G to every other node in G, i.e. create a complete graph. ----- HASTA AQUI adj_list1.cpp
3) Create a new node i.
4) Pick a node j uniformly at random from the graph G. Set P = (k(j)/k_tot)^a.
5) Pick a real number R uniformly at random between 0 and 1.
6) If P > R then add j to i's adjacency list.
7) Repeat steps 4 - 6 until i has m nodes in its adjacency list.
8) Add i to the adjacency list of each node in its adjacency list.
9) Add i to to the graph.
10) Repeat steps 3 - 9 until there are N nodes in the graph.
*/

int main(){

	Graph g;
	srand(time(0));

	ifstream fileIn; //input-stream variable
	//ofstream output_file;
	
	
	int n_S=0,n_I=0,n_D=0;//number of total S,I,D.
	int status_for_switch, node_identifier;
	float angle, t=0,dt,flag, vel;
	int tf, radius, n, BoxX, BoxY, num_sim;
	float beta;
	int ni = 4, nf = 100;
	//Person P[2*n], tmp[2*n];

	//cout << "Write the name of the file you want to open: ";
	//cin >> inputfilename;

	fileIn.open("datin.txt");
	if(!fileIn)
	{
		cout << "Couldn't find the input file" << endl;
		exit(-1);
	}

	fileIn >> num_sim >> beta >> n >> BoxX >> BoxY >> tf >> dt >> radius;

	fileIn.close();

	cout << "simulaciones: " << num_sim << ", "
		<< "beta: " << beta << ", "
		<< "N: " << n << ", "
		<< "BoxX: " << BoxX << ", "
		<< "BoxY: " << BoxY << ", "
		<< "tf: " << tf << ", "
		<< "dt: " << dt << ", "
		<< "radius: " << radius << ", " << endl;

	tostr(beta);

	int simulations = 1;
	int sim_number;

	for (sim_number = 0; sim_number < simulations; sim_number++)
	{	
		//cout < sim_number < endl;
		
		g.clearGraph();

	// Initial nodes
		for (int i = 0; i < ni; ++i)
		{
			Node v;
			v.setID(i);
			g.addNode(v);
		}

	// Now we need to connect them somehow.
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < ni; ++j)
			{
				if (i != j)
				{
					g.addEdgeByID(i,j); //(from, to)
				}
			}

			//cout << "NumberOfEdgesByNode: " << g.NumberOfEdgesByNode(i) << endl;

		}

		//g.printGraph();


	// Add a new node
		int contador = 0;
		int new_nodes = ni;
		int sum_ki, kj = 3; // ki -> degree of old nodes... kj -> degree of new node.
		int random_node=0;

		while(contador != nf)
		{
			// Add a new node
			Node v;
			v.setID(new_nodes);
			g.addNode(v);
			//---------------

			// busca la probabilidad de cada nodo existente y los 3 mayores se acoplan.
			// I am wondering: does a power law distribution already implies that there are hubs in the network?
			// Degree k is the number of Edges P(i,j) = kj/SUM(ki)
			//Hacer funcion que calcule la probabilidad de cada nodo!
			
			// Sum the edges of all existing nodes.
			sum_ki = 0;
			for (int i = 0; i < new_nodes; ++i)
			{
				sum_ki += g.NumberOfEdgesByNode(i); 
			}

			//set <int> visitados;
			double Pk = 0;
			float dice = 0;
			int conections_per_new_node = 0;
			// Pick a random node, check the prob to connect and throw a random number, if its greater than the prob, connect it.
			/*do{
				random_node = rand() % new_nodes;
				Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1.22); //why 1.22?
				cout << "Pk: " << Pk << endl;
				dice = (float)rand()/RAND_MAX;
				//cout << dice << endl;
				if (Pk > dice)
				{
					g.addEdgeByID(new_nodes,random_node);
					conections_per_new_node++;
				}
					
			}while(conections_per_new_node < kj);
			*/
			while(conections_per_new_node < kj)
			{
				random_node = rand() % new_nodes;
				Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1.22); //why 1.22?
				cout << "Pk: " << Pk << endl;
				dice = (float)rand()/RAND_MAX;
				//cout << dice << endl;
				if (Pk > dice && g.checkIfEdgeExistsByID(new_nodes,random_node)==false)
				{
					g.addEdgeByID(new_nodes,random_node);
					conections_per_new_node++;
				}
			}
			contador++;
			new_nodes++;
		}
		g.printGraph();

		// This section is designed to plot the characteristic curve of the SF network.
		int ListNumberOfEdges[nf+1];
		int FreqDist[nf+1];

		for (int i = 0; i < nf+2; ++i)
		{
			FreqDist[i] = 0;
		}

		for (int i = 0; i < nf+1; ++i)
		{
			ListNumberOfEdges[i] = g.NumberOfEdgesByNode(i); // Node "i" has "g.NumberOfEdgesByNode(i)" Edges.
			//cout << ListNumberOfEdges[i] << endl;
			FreqDist[ListNumberOfEdges[i]]++; //Now ListNumberOfEdges[i] is the index of the FreqDist array and we add one to that index every time we find it.
			//cout << FreqDist[ListNumberOfEdges[i]] << endl;
		}
		// Printing section
		//--------------------------------
		std::string output = "sf2_network_" + std::to_string(sim_number);
		std::ofstream sf(output);

		//output_file.open("sf_network");

		for (int i = 0; i < nf+1; ++i)
		{
			if (i >= 0)
			{

				sf << i << " " << FreqDist[i] << endl;
			}
		}

		//output_file.close();
		//--------------------------------
		sf.close();
	}
	return 0;
}//end

void tostr(float beta){
	betastr = int(beta*10);
}