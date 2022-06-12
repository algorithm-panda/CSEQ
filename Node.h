/* 
 * File:   Node.h
 * Author: hjf
 *
 * Created on April 15, 2016, 9:11 PM
 */

#ifndef NODE_H
#define	NODE_H
#include<iostream>
#include<string>
#include<vector>
#include<set>
using namespace std;

//#define alpha 0.5
#define BASELINE 0
#define ONLY_ONE_EX 1
#define Error_bound 1e-9

class Node{
public:
    int business_id;
    string name;
    double latitude;
    double longitude;
    double stars;
    double review_count;
    double type_arr[20];
    set<string> categories;
};

class StringNode{
public:
	string business_id;
    string name;
    string latitude;
    string longitude;
    string stars;
    string review_count;
    vector<string> categories;
};

class PNode{
public:    
    int id;
    double weight;
    double distance;
    PNode(){}
    PNode(int i, double w, double dist)
    {
        id=i;
        weight=w;
        distance=dist;
    }
};

bool comp(PNode &A, PNode& B)
{
    return A.weight > B.weight;
}

bool comp2(PNode &A, PNode& B){
	return A.distance < B.distance;
}

class QueueNode{
public:
		vector<int> id;
		double weight;
                double avgDist;
		QueueNode();
               
                
		QueueNode(vector<int> &id, double weight, double dist)
                {
                    this->id=id;
                    this->weight=weight;
                    avgDist=dist;
                }
		friend bool operator < (QueueNode a, QueueNode b)  //递增for priority_queue,   递减for set
		{
			return a.weight > b.weight || (a.weight == b.weight && a.avgDist< b.avgDist);
		}
                
};

#endif	/* NODE_H */

