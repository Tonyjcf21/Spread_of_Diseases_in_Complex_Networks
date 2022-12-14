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

using namespace std;

#define TAU 6.283185307

// g++ sf_network_rewire_SIR.cpp -o sfn -std=c++2a

class Edge;
class Node;

class Edge
{
	public:
		int destinationID;
		int n_edges = 0;

		Edge() {}
		Edge(int destID){
			destinationID = destID;
		}

		//void setEdgeValues(int destID){
		//	destinationID = destID;
		//}

		int getDestinationID(){
			return destinationID;
		}

};

class Node {

public:
	int node_id;
	int node_status; //0=S, 1=I, 2=D
	int recovery_time;
	int dead_time;
	list <Edge> edgeList;

	Node(){
		node_id = 0;
		node_status = 0;
		recovery_time = 10;
		dead_time = 2;
	}

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

	void setStatus(int status){
		node_status = status;
	}

	int getStatus(){
		return node_status;
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

    void clearEdges(){
    	edgeList.clear();
    }
};

class Graph{

	vector<Node> nodes;

public:

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
	    }
    }

    void PrepareToRewire(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	temp.clearEdges();
    }

    void NodeInfection(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	temp.setStatus(1);
    }

    void NodeRecovery(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	temp.setStatus(0);
    }

    int* EdgesOfNode(int node){
    	Node v = getNodeByID(node);
    	list < Edge > e;
    	e = v.getEdgeList();
    	int* arr = new int[2000];
    	int i = 0;

    	for (auto it = e.begin(); it != e.end(); it++) {
    		arr[i] = it -> getDestinationID();
    		i++;
    	}

    	return arr;
    }

/*
    int NodeNumber(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	return temp;
    }
    */
};


int main(){

	Graph g;
	srand(time(0));

	int simulations = 1;
	int sim_number;
	int m = 5; // ki -> degree of old nodes... m -> degree of new node.
	int m0 = m+1, N = 100;

	for (sim_number = 0; sim_number < simulations; sim_number++)
	{	
		g.clearGraph();

		int contador = 0;
		int sum_ki; 
		int new_nodes = m0;
		int random_node;
		double Pk = 0;
		float dice = 0;

		// Initial nodes
		for (int i = 0; i < m0; ++i)
		{
			Node v;
			v.setID(i);
			g.addNode(v);
		}
	
		// Now we need to connect them somehow. This is fully connected.
		for (int i = 0; i < m0; ++i)
		{
			for (int j = 0; j < m0; ++j)
			{
				if (i != j)
				{
					g.addEdgeByID(i,j); //(from, to)
				}
			}
			//cout << "NumberOfEdgesByNode: " << g.NumberOfEdgesByNode(i) << endl;
		}
		g.printGraph();

		// Add a new node
		while(contador != N)
		{
			cout << contador << " ";
			// Add a new node
			Node v;
			v.setID(new_nodes);
			g.addNode(v);
			//---------------

			// I am wondering: does a power law distribution already implies that there are hubs in the network?
			
			// Sum the edges of all existing nodes.
			sum_ki = 0;
			for (int i = 0; i < new_nodes; ++i)
			{
				sum_ki += g.NumberOfEdgesByNode(i); //sum_ki is twice the total edges.
			}

			int conections_per_new_node = 0;

			// Pick a random node, check the prob to connect and throw a random number, if its lower than the prob, connect it.
			while(conections_per_new_node < m)
			{
				random_node = rand() % new_nodes;
				Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1); //why 1.22? Prob that the new node connects to selected node.
				
				dice = (float)rand()/RAND_MAX;
				
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

		//----------------------------------------------------------------------------------------------------------------
		//                                             Rewiring section
		//----------------------------------------------------------------------------------------------------------------
		int amount_of_nodes_to_rewire = 5;
		float choice;
		float rewiring_rate[N+1];
		int rewiring_node, numberOfEdges, contador2 = 0;

		for (int i = 0; i < N; ++i)
		{
			rewiring_rate[i] = m/g.DegreeOfNodeToAttach(i);
		}

		while(contador2 < amount_of_nodes_to_rewire)
		{
			bool flag = false;

			while(flag != true)
			{
				rewiring_node = rand() % N;
				choice = (float)rand()/RAND_MAX;

				if(choice < rewiring_rate[rewiring_node])
				{
					numberOfEdges = g.NumberOfEdgesByNode(rewiring_node);
					g.PrepareToRewire(rewiring_node);
					flag = true;
					contador2++;
					cout << "node " << rewiring_node << " has " << numberOfEdges << " edges." << endl;
				}
			}

			int conections_per_rewired_node = 0;

			// Pick a random node, check the prob to connect and throw a random number, if its lower than the prob, connect it.
			while(conections_per_rewired_node < numberOfEdges)
			{
				random_node = (rand() % (N-1));
				if (rewiring_node != random_node)
				{
					Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1); 

					dice = (float)rand()/RAND_MAX;
					
					if (Pk > dice && g.checkIfEdgeExistsByID(rewiring_node,random_node)==false)
					{
						g.addEdgeByID(rewiring_node,random_node);
						conections_per_rewired_node++;
					}
				}
			}
			cout << "node " << rewiring_node << " has been rewired!" << endl;
		}

		//----------------------------------------------------------------------------------------------------------------
		//                  This section is designed to plot the characteristic curve of the SF network.
		//----------------------------------------------------------------------------------------------------------------
		int ListNumberOfEdges[N+1];
		int FreqDist[N+1];

		// We need to initialize the array to avoid getting random elements from accupied memory spaces.
		for (int i = 0; i < N+2; ++i)
		{
			FreqDist[i] = 0;
		}

		for (int i = 0; i < N+1; ++i)
		{
			ListNumberOfEdges[i] = g.NumberOfEdgesByNode(i); // Node "i" has "g.NumberOfEdgesByNode(i)" Edges.
			FreqDist[ListNumberOfEdges[i]]++; //Now ListNumberOfEdges[i] is the index of the FreqDist array and we add one to that index every time we find it.
		}

		//----------------------------------------------------------------------------------------------------------------
		//                                             Printing section
		//----------------------------------------------------------------------------------------------------------------
		std::string output = "sf2_network10000_rewired500_" + std::to_string(sim_number);
		std::ofstream sf(output);

		for (int i = 0; i < N+1; ++i)
		{
			if (i >= m)
			{
				sf << i << " " << FreqDist[i] << endl;
			}
		}
	
		//----------------------------------------------------------------------------------------------------------------
		
		sf.close();
	}

	int t = 0, tf = 1;
	float odds, beta = 0.6;
	list<int> infecteds;


	// Initializing the disease.
	int random_infection = rand() % N;
	g.NodeInfection(random_infection);

	// Finding the edges of the infected node.
	int* Connections = g.EdgesOfNode(random_infection);

	// Infecting the neighbours of the infected node with a probability.
	for (int i = 0; i < g.NumberOfEdgesByNode(random_infection); ++i)
	{
		odds = (float)rand()/RAND_MAX;
		if (odds < beta)
		{
			g.NodeInfection(Connections[i]);
			infecteds.push_back(Connections[i]);
		}
	}

	// The initial infection is recovered.
	g.NodeRecovery(random_infection);

	// Remove data from the dynamic array.
	delete[] Connections;

	//Now that we have a list of infecteds (Hopefully), we get the loops running.
	while(t<tf)
	{
		int count = 0;

		for (auto it = infecteds.begin(); it != infecteds.end(); it++)
		{
			// This will be the infected node that is spreading the disease among its edges.
			int specific_node = *(it); // list iterator to int conversion.
			int* Connections = g.EdgesOfNode(specific_node); // Saving the edges from this node (specific_node).

			cout << "trying a print here..." << endl;
			cout << specific_node << " " << Connections[0] << " " <<Connections[1] << endl;

			// Infecting nodes connected to specific_node with a certain probability.
			for (int j = 0; j < g.NumberOfEdgesByNode(Connections[count]); ++j)
			{
				odds = (float)rand()/RAND_MAX;
				if (odds < beta)
				{
					g.NodeInfection(Connections[j]); // Infecting the node.
					infecteds.push_back(Connections[j]); // Adding the node to the infected list.
				}
			}
			count++;
			g.NodeRecovery(specific_node); // Now that this node has infected some of its neighbours. It has recovered.
			delete[] Connections; // Cleaning the edge list of this node so we can put the list of the next node.
		}
		t++;
	}

	return 0;
}//end