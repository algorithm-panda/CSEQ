
#include <cstdlib>
#include <cstdio>
#include<fstream>
#include<vector>
#include<set>
#include<sstream>
#include<cmath>
#include <queue>
#include<assert.h>
#include<time.h>
#include <algorithm>
#include <sys/time.h>
#include<map>
#include<cstring>
#include<algorithm>
#include <unordered_map>
#include "Node.h"
#include "utility.h"
#include <iostream>


using namespace std;
int sumHash[1000000];
vector<vector<int> > detailHash;
vector<vector<int> > changedDetailHash;
vector<Node> NodeSet;
vector<double> RadioSet;
int numMapping[10];
int COUNT = 0;

double percentage = 0;
double total_possible = 0;
vector<vector<double> > SpatialVector;

map<string, set<int> > InvertedList;

double alpha = 0.5;
double magnitude;
int K = 5;
int LatSplitNum;
int LonSplitNum;
double latLength;
double lonLength;
double rate_p = 1.5;
double Error_bound_1 = 1e-9;

double dfs_time = 0;

double catePairDistance[10][10];

unordered_map<int, double> itemPairDistance[10][10];
unordered_map<int, double> searchPairDistance[10][10];
unordered_map<int, int> searchPairDistanceIndex[10][10];
unordered_map<int,double> prefixItemPairDistance[10][10];

double CostInSpatialUb = 0;
double CostInPreSort = 0;
double CostInSimMap = 0;
double CostInPreNear = 0;

vector<double> minLat;
vector<double> minLon;
vector<double> maxLat;
vector<double> maxLon;



unordered_map<int,int> gridDict[51][15];
vector< unordered_map<int,vector<PNode> > > gridListDict[51];


bool MoveFlag;

long double sim(Node a, Node b);
long double avgTextualSim(Node a, int index, vector<vector<int> >&QuerySet);
bool checkDist(double R, int a, int b);
vector<PNode> FilterTypeByText(vector<Node> &nodeset, string type, double,
                               double);
vector<string> split(string s, char c) {
    vector<string> ans;
    string temp = "";
    for (int i = 0; i < s.size(); i++) {
        if (s[i] == c) {
            if (temp.size() > 0 && temp[0] != ',')
                ans.push_back(temp);
            else if (temp.size() == 0)
                ans.push_back(temp);
            temp = "";
        } else if (s[i] != ':' && s[i] != ' ' && s[i] != ',') {
            temp += s[i];
        }
    }
    if (temp.size() > 0)
        ans.push_back(temp);
    return ans;
}

void stringToint(int &a, string s) {
    std::stringstream ss;
    ss << s;
    ss >> a;
}

void stringToDouble(double &a, string s) {
    std::stringstream ss;
    ss << s;
    ss >> a;
}


void BuildIndex() {

    int num = 0;
    int n = NodeSet.size();
    for (int i = 0; i < n; i++) {
        for (string e : NodeSet[i].categories) {
            if (InvertedList.find(e) != InvertedList.end()) {
                InvertedList[e].insert(i);
            } else {
                set<int> s;
                s.insert(i);
                InvertedList.insert(pair<string, set<int> >(e, s));
            }
        }
    }
}

void Read(vector<Node> &nodeset, string path) {
    ifstream ifile;
    ifile.open(path.c_str());
    cout << path.c_str() << endl;
    if (ifile.fail()) {
        cout << "file cannot found" << endl;
    }
    string str;

    double maxstar = 0;
    double maxReviewCount = 0;
    int tempCnt = 0;
    int line = 0;

    while (getline(ifile,str)) {
        line++;
        Node node;
        vector<string> ans = split(str, ',');

        node.business_id = line;
        node.name = ans[0].substr(0,10);
        stringToDouble(node.latitude, ans[1]);
        stringToDouble(node.longitude, ans[2]);
        for (int j = 0;j < 20; j++) {
            stringToDouble(node.type_arr[j], ans[j+4]);
        }


        node.categories.insert(ans[3]);

        nodeset.push_back(node);
    }

    ifile.close();


}
double dist(double lat1, double lon1, double lat2, double lon2) {
    return Spherical_distance(lat1, lon1, lat2, lon2);
}

double dist(vector<Node>& nodeset, int id1, int id2) {
    return Spherical_distance(nodeset[id1].latitude, nodeset[id1].longitude,
                              nodeset[id2].latitude, nodeset[id2].longitude);
}

vector<PNode> FilterType(vector<Node> &nodeset, int id, int index,
                         vector<vector<int> > &querySet, double r) {
    vector<PNode> ans;

    int count = 0;


    bool used[nodeset.size() + 1];
    memset(used, 0, sizeof(used));
    for (auto e : nodeset[id].categories) {
        for (auto item : InvertedList[e]) {
            if (used[item] || checkDist(r, id, item) == false)
                continue;
            ans.push_back(
                    PNode(item, avgTextualSim(nodeset[item], index, querySet),
                          dist(NodeSet[id].latitude, NodeSet[id].longitude,
                               NodeSet[item].latitude,
                               NodeSet[item].longitude)));
            used[item] = 1;

        }
    }
    return ans;
}


vector<PNode> FilterTypeByText(vector<Node> &nodeset, string type,
                               double review_count, double stars) {
    double nearZero = 0.000001;
    vector<PNode> ans;
    for (int item : InvertedList[type]) {
        if (review_count > nearZero) {
            if (nodeset[item].review_count < review_count)
                continue;
        }

        if (stars > nearZero) {
            if (nodeset[item].stars < stars)
                continue;
        }
        ans.push_back(PNode(item, 0, 0));
    }

    return ans;
}







bool check(vector<Node> &NodeSet, int id1, double dist1, int id2, double dist2) //j dominates i
{

    if (dist2 > dist1)
        return false;

    if (NodeSet[id2].stars < NodeSet[id1].stars)
        return false;

    if (NodeSet[id2].review_count < NodeSet[id1].review_count)
        return false;

    return true;
}




long double sim(Node a, Node b) {
    long double ans = 0;
    double fenzi = 0;
    double fenmu_1 = 0;
    double fenmu_2 = 0;
    for (int i = 0;i < 20; i++) {
        fenzi += a.type_arr[i] * b.type_arr[i];
        fenmu_1 += a.type_arr[i] * a.type_arr[i];
        fenmu_2 += b.type_arr[i] * b.type_arr[i];
    }

    if (fenmu_1 == 0 || fenmu_2 == 0)
        return 0;

    return fenzi / (sqrt(fenmu_1) * sqrt(fenmu_2));

}

long double avgTextualSim(Node a, int index, vector<vector<int> >&QuerySet) {
    int n = QuerySet.size();
    long double ans = 0;
    for (int i = 0; i < n; i++) {
        long double x = sim(a, NodeSet[QuerySet[i][index]]);


        ans += x;
    }

    return ans / n;
}



long double GetTextualSim(vector<int> & Query, vector<int>& combination) {
    long double ans = 0;
    int n = Query.size();

    for (int i = 0; i < n; i++) {
        ans += sim(NodeSet[Query[i]], NodeSet[combination[i]]);
    }

    return ans / n;
}

vector<double> getSpatialVector(vector<int> & y) {
    vector<double> ans;
    int n = y.size();
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++) {
            ans.push_back(
                    dist(NodeSet[y[i]].latitude, NodeSet[y[i]].longitude,
                         NodeSet[y[j]].latitude, NodeSet[y[j]].longitude));
        }
    return ans;
}

long double GetSpatialSim(vector<int> & Query, vector<int>& combination) {
    long double ans = 0;
    vector<double> V1;
    vector<double> V2;

    V1 = getSpatialVector(Query);
    V2 = getSpatialVector(combination);

    for (int i = 0; i < V1.size(); i++)
        ans += V1[i] * V2[i];

    long double X = 0, Y = 0;
    for (int i = 0; i < V1.size(); i++) {
        X += (V1[i] * V1[i]);
        Y += (V2[i] * V2[i]);
    }
    ans /= (sqrt(X) * sqrt(Y));

    return ans;
}

long double GetSpatialSim(int index, vector<int>& combination) {
    long double ans = 0;
    vector<double> V2;

    V2 = getSpatialVector(combination);

    for (int i = 0; i < SpatialVector[index].size(); i++)
        ans += SpatialVector[index][i] * (V2[i]);

    long double X = 0, Y = 0;
    for (int i = 0; i < SpatialVector[index].size(); i++) {
        X += (SpatialVector[index][i] * SpatialVector[index][i]);
        Y += (V2[i] * V2[i]);
    }
    ans /= sqrt(Y);
    ans /= sqrt(X);

    return ans;
}


double getDimVector(vector<int> & comb, int kcnt) {
    double ans = 0;
    for (int i = 0;i < comb.size(); i++) {
        int nodeId = comb[i];
        double latLength = (maxLat[i] - minLat[i]) / kcnt;
        double lonLength = (maxLon[i] - minLon[i]) / kcnt;
        int x = int((NodeSet[nodeId].latitude - minLat[i]) / latLength);
        int y = int((NodeSet[nodeId].longitude - minLon[i]) / lonLength);
        ans = max(ans, dist(x * latLength + minLat[i], y * lonLength + minLon[i], (x+1) * latLength + minLat[i], (y+1) * lonLength + minLon[i]));
    }
    return ans;
}

double get_approx_bound(vector<int> & comb, vector<Node> &NodeSet, int kcnt) {
    for (int i = 0;i < comb.size(); i++) {
        int nodeId = comb[i];
        double latLength = (maxLat[i] - minLat[i]) / kcnt;
        double lonLength = (maxLon[i] - minLon[i]) / kcnt;
        int x = int((NodeSet[nodeId].latitude - minLat[i]) / latLength);
        int y = int((NodeSet[nodeId].longitude - minLon[i]) / lonLength);

    }

}

long long calGridCnt(vector<int> & comb, int kcnt) {
    long long ans = 1;
    for (int i = 0;i < comb.size(); i++) {
        int nodeId = comb[i];
        double latLength = (maxLat[i] - minLat[i]) / kcnt;
        double lonLength = (maxLon[i] - minLon[i]) / kcnt;
        int x = int((NodeSet[nodeId].latitude - minLat[i]) / latLength);
        int y = int((NodeSet[nodeId].longitude - minLon[i]) / lonLength);
        int gridId = x * kcnt + y;
        ans *= gridDict[i][kcnt][gridId];
    }
    return ans;
}


double GetSpatialUPSim(int index, vector<int>& combination, int level) {
    double ans = 0;
    vector<double> V2;

    V2 = getSpatialVector(combination);
    double dim = getDimVector(combination, level);


    for (int i = 0; i < SpatialVector[index].size(); i++)
        ans += SpatialVector[index][i] * (V2[i] + dim);

    long double X = 0, Y = 0;
    for (int i = 0; i < SpatialVector[index].size(); i++) {
        X += (SpatialVector[index][i] * SpatialVector[index][i]);
        Y += (max(0.0,V2[i] - dim) * max(0.0, V2[i] - dim));
    }
    ans /= sqrt(Y);
    ans /= sqrt(X);

    return min(1.0,ans);
}


long double sim(vector<int>& combination1, vector<int>& combination2) {
    return alpha * GetTextualSim(combination1, combination2)
           + (1 - alpha) * GetSpatialSim(combination1, combination2);
}

long double spatialSim(vector<vector<int> > &QuerySet, vector<int>& combination,
                       double Currentweight) {
    long double ans = 0;

    for (int i = 0; i < QuerySet.size(); i++) {
        ans += alpha * GetSpatialSim(i, combination);
    }
    return  ans / QuerySet.size();

}

long double distanceRate(vector<vector<int> > &QuerySet, vector<int>& combination,
                         double Currentweight) {
    vector<double> V1,V2;
    V1.clear();
    V2.clear();

    for (int i = 0;i < QuerySet.size(); i++) {
        V1 = getSpatialVector(QuerySet[i]);
        V2 = getSpatialVector(combination);

        long double X = 0, Y = 0;
        for (int j = 0; j < V1.size(); j++) {
            X += (V1[j] * V1[j]);
            Y += (V2[j] * V2[j]);
        }
        X = sqrt(X);
        Y = sqrt(Y);
        return X/Y;
    }

}


double UPsim(vector<vector<int> > &QuerySet, vector<int>& combination,
             double Currentweight, int level) {
    double ans = 0;

    for (int i = 0; i < QuerySet.size(); i++) {
        ans += alpha * GetSpatialUPSim(i, combination, level);
    }
    return (1 - alpha) * Currentweight / QuerySet[0].size()
           + ans / QuerySet.size();

}

long double sim(vector<vector<int> > &QuerySet, vector<int>& combination,
                double Currentweight) {
    long double ans = 0;

    for (int i = 0; i < QuerySet.size(); i++) {
        ans += alpha * GetSpatialSim(i, combination);
    }
    return (1 - alpha) * Currentweight / QuerySet[0].size()
           + ans / QuerySet.size();

}

long double sim(vector<vector<int> > &QuerySet, vector<int>& combination) {
    long double ans = 0;

    for (int i = 0; i < QuerySet.size(); i++) {
        ans += (1 - alpha) * GetTextualSim(QuerySet[i], combination)
               + alpha * GetSpatialSim(QuerySet[i], combination);
    }
    return ans * 1.0 / QuerySet.size();

}

double GetAvgDist(vector<Node> &NodeSet, vector<int> &combination) {
    double ans = 0;
    int n = combination.size();
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            ans += dist(NodeSet[combination[i]].latitude,
                        NodeSet[combination[i]].longitude,
                        NodeSet[combination[j]].latitude,
                        NodeSet[combination[j]].longitude);

    return 2 * ans / (n * (n - 1));
}

int getGridId(int nodeId,double minLat, double minLon, double maxLat, double maxLon, int cnt) {
    int dcnt = cnt;
    double latLength = (maxLat - minLat + 1e-9) / dcnt;
    double lonLength = (maxLon - minLon + 1e-9) / dcnt;
    int x = int((NodeSet[nodeId].latitude - minLat) / latLength);
    int y = int((NodeSet[nodeId].longitude - minLon) / lonLength);


    return x * dcnt + y;
}





unordered_map<int,pair<double, double> > gridSpatialSimMap[10];



unordered_map<int,vector<PNode> > gridList[10];
unordered_map<int,vector<PNode> > gridspaList[10];


void dfs(vector<unordered_map<int,vector<PNode> > > &List, int num, vector<vector<int> > &querySet,
         vector<int> & combination, priority_queue<QueueNode> & que,
         double CurrentWeight,double totminLat, double totmaxLat, double totminLon, double totmaxLon,
         vector<double> & region_max_value) {
    if (num == querySet[0].size()) {
        priority_queue<pair<double,vector<int> > > q;
        double sum = 0;
        set<vector<int> > ps;
        ps.clear();
        vector<int> sum_vec;
        sum_vec.clear();
        for (int i = 0;i < num; i++) {

            sum += List[i][combination[i]][0].weight;
            sum_vec.push_back(0);
        }
        q.push(make_pair(sum, sum_vec));

        for (int i = 0;i < K; i++) {

            if (q.empty()) break;
            pair<double,vector<int> > wc = q.top();
            q.pop();

            COUNT++;

            vector<int> candi_vector;
            candi_vector.clear();
            bool in_flag = false;
            for (int j = 0;j < num; j++) {
                candi_vector.push_back(List[j][combination[j]][wc.second[j]].id);
                if (j == 0) {
                    double lat1 = NodeSet[List[j][combination[j]][wc.second[j]].id].latitude;
                    double lon1 = NodeSet[List[j][combination[j]][wc.second[j]].id].longitude;
                    if (lat1 >= totminLat && lat1 < totmaxLat && lon1 >= totminLon && lon1 < totmaxLon) {
                        in_flag = true;
                    }
                }
                if (j > 0 && List[j][combination[j]][wc.second[j]].id == 0) {

                }

            }
            double similarity = sim(querySet, candi_vector, wc.first);


            double rate = distanceRate(querySet, candi_vector, wc.first);



            if (que.size() >= K && ((1 - alpha) * wc.first / querySet[0].size() + alpha * 1.0) < que.top().weight + Error_bound_1)
                return;

            if ((rate > rate_p) or (rate < 1.0 / rate_p) or (not in_flag)) {
                if (rate > ((1 + rate_p) * rate_p)) {
                    return;
                }

                i -= 1;
            } else {

                if (que.size() < K) {
                    que.push(
                            QueueNode(candi_vector, similarity,
                                      GetAvgDist(NodeSet, candi_vector)));

                } else if (similarity > que.top().weight) {
                    que.pop();
                    que.push(
                            QueueNode(candi_vector, similarity,
                                      GetAvgDist(NodeSet, candi_vector)));


                }
            }
            int num1 = 0;
            if (not in_flag) {
                num1 = 1;
            } else  num1 = num;
            for (int j = 0; j < num1; j++) {
                if (wc.second[j] + 1 < List[j][combination[j]].size()) {
                    double sum_temp = wc.first - List[j][combination[j]][wc.second[j]].weight +
                                      List[j][combination[j]][wc.second[j] + 1].weight;
                    wc.second[j]++;
                    if (ps.find(wc.second) == ps.end()) {
                        q.push(make_pair(sum_temp, wc.second));
                        ps.insert(wc.second);
                    }
                    wc.second[j]--;
                }
            }

        }

        return;
    }

    int cSize = List.size();

    for (auto &gridId : List[num]) {
        combination[num] = gridId.first;

        if (que.size() >= K) {
            double UB = 0;
            double cw = 0;
            for (int j = 0;j <= num; j++) {
                cw += List[j][combination[j]][0].weight;
            }
            for (int j = num+1; j < querySet[0].size(); j++) {
                cw += region_max_value[j];
            }
            UB = (1-alpha) * cw / querySet[0].size() + alpha * 1.0;
            if (UB < min(Error_bound_1 + que.top().weight,1.0)) {

                continue;
            }
        }

        dfs(List, num + 1, querySet, combination, que,
            CurrentWeight + 0,totminLat, totmaxLat, totminLon, totmaxLon, region_max_value);
    }
}





void splitDFS(vector<Node> &NodeSet, vector< vector<PNode> > &NodeListSet, vector<vector<int> > &querySet, double length_limit, int D, int odd, vector<int> & combination, priority_queue<QueueNode> & que, long long a1, long long a2,
              double totminLat, double totmaxLat, double totminLon, double totmaxLon) {

    double newminLat = 1e18;
    double newminLon = 1e18;
    double newmaxLat = -1e18;
    double newmaxLon = -1e18;

    double minLatList[10];
    double minLonList[10];
    double maxLatList[10];
    double maxLonList[10];

    int minmaxI = 0;
    for (auto nodeList:NodeListSet) {
        minLatList[minmaxI] = 1e18;
        minLonList[minmaxI] = 1e18;
        maxLatList[minmaxI] = -1e18;
        maxLonList[minmaxI] = -1e18;
        for (auto p:nodeList) {
            newminLat = min(newminLat, NodeSet[p.id].latitude);
            newminLon = min(newminLon, NodeSet[p.id].longitude);
            newmaxLat = max(newmaxLat, NodeSet[p.id].latitude);
            newmaxLon = max(newmaxLon, NodeSet[p.id].longitude);
            minLatList[minmaxI] = min(minLatList[minmaxI], NodeSet[p.id].latitude);
            minLonList[minmaxI] = min(minLonList[minmaxI], NodeSet[p.id].longitude);
            maxLatList[minmaxI] = max(maxLatList[minmaxI], NodeSet[p.id].latitude);
            maxLonList[minmaxI] = max(maxLonList[minmaxI], NodeSet[p.id].longitude);
        }
        minmaxI += 1;
    }
    double midLat = (totminLat + totmaxLat) / 2.0;
    double midLon = (totminLon + totmaxLon) / 2.0;
    double distRan1 = dist(midLat, midLon, midLat, totmaxLon);
    double distRan2 = dist(midLat, midLon, totmaxLat, midLon);




    if ((distRan1 <= length_limit * rate_p + Error_bound_1) or (distRan2 <= length_limit * rate_p + Error_bound_1)) {
        vector< unordered_map<int,vector<PNode> > > gridList;
        gridList.clear();
        gridList.resize(querySet[0].size());
        double totdis = dist(newminLat,newminLon, newmaxLat, newmaxLon);

        double rr = min(2.0,(2*sqrt(3)* rate_p * totdis / D / length_limit + 1));

        int type_cnt = 0;
        vector<double> region_max_value;
        region_max_value.clear();
        for (auto nodeList: NodeListSet) {
            double max_att_val = 0;
            for (auto p:nodeList) {
                int gridId = getGridId(p.id, minLatList[type_cnt], minLonList[type_cnt], maxLatList[type_cnt], maxLonList[type_cnt], D);
                if (gridList[type_cnt].find(gridId) == gridList[type_cnt].end()) {
                    gridList[type_cnt][gridId].clear();
                }
                if (gridList[type_cnt][gridId].size() < K * rr) {
                    gridList[type_cnt][gridId].push_back(p);
                }
                max_att_val = max(max_att_val, p.weight);
            }
            type_cnt += 1;
            region_max_value.push_back(max_att_val);

        }

        clock_t s_dfs = clock();
        dfs(gridList, 0, querySet, combination, que, 0, totminLat, totmaxLat, totminLon, totmaxLon, region_max_value);
        clock_t e_dfs = clock();
        dfs_time += (e_dfs - s_dfs ) /  (double) CLOCKS_PER_SEC;

    }
    else {
        vector<vector<PNode> > leftListSet;
        vector<vector<PNode> > rightListSet;
        leftListSet.clear();
        rightListSet.clear();

        bool leftFlag = true;
        bool rightFlag = true;
        for (auto pList:NodeListSet) {
            vector<PNode> leftList;
            vector<PNode> rightList;
            leftList.clear();
            rightList.clear();
            for (auto p:pList) {
                if (odd % 2 == 0) {
                    double dis1 = dist(midLat, NodeSet[p.id].longitude, NodeSet[p.id].latitude, NodeSet[p.id].longitude) ;


                    if ((NodeSet[p.id].latitude < midLat) or (dis1 < length_limit * rate_p)) {
                        leftList.push_back(p);
                    }

                    if ((NodeSet[p.id].latitude >= midLat) or (dis1 < length_limit * rate_p)) {
                        rightList.push_back(p);
                    }

                } else {

                    double dis1 = dist(NodeSet[p.id].latitude, midLon, NodeSet[p.id].latitude, NodeSet[p.id].longitude) ;


                    if ((NodeSet[p.id].longitude < midLon) or (dis1 < length_limit * rate_p)) {
                        leftList.push_back(p);
                    }

                    if ((NodeSet[p.id].longitude >= midLon) or (dis1 < length_limit * rate_p)) {
                        rightList.push_back(p);
                    }
                }
            }

            if (leftList.size() == 0) {
                leftFlag = false;
            }
            if (rightList.size() == 0){
                rightFlag = false;
            }
            leftListSet.push_back(leftList);
            rightListSet.push_back(rightList);
        }

        if (leftFlag) {
            if (odd % 2 == 0) {
                splitDFS(NodeSet, leftListSet, querySet, length_limit, D, (odd + 1), combination, que, a1 * 2 + 0,
                         a2 * 2, totminLat, midLat, totminLon, totmaxLon);
            } else {
                splitDFS(NodeSet, leftListSet, querySet, length_limit, D, (odd + 1), combination, que, a1 * 2 + 0,
                         a2 * 2, totminLat, totmaxLat, totminLon, midLon);
            }
        }
        if (rightFlag) {
            if (odd % 2 == 0) {
                splitDFS(NodeSet, rightListSet, querySet, length_limit, D, (odd + 1), combination, que, a1 * 2 + 1,
                         a2 * 2, midLat, totmaxLat, totminLon, totmaxLon);
            } else {
                splitDFS(NodeSet, rightListSet, querySet, length_limit, D, (odd + 1), combination, que, a1 * 2 + 1,
                         a2 * 2, totminLat, totmaxLat, midLon, totmaxLon);
            }
        }
    }
}



vector<QueueNode> getTopK(vector<Node> &NodeSet, vector<vector<int> > &querySet,
                          double r,int D, int K) {
    priority_queue<QueueNode> que;
    priority_queue<QueueNode> UPque;
    vector<vector<PNode> > oriList;
    vector<vector<PNode> > oriListSpatial;
    vector<vector<PNode> > List[50];
    vector<vector<PNode> > upList[50];
    vector<vector<PNode> > ListSpatial[50];
    vector<vector<PNode> > upListSpatial[50];
    vector<vector<int> > changedQuerySet;

    long long possible = 1;
    double data_ratio = 1;
    for (int i = 0; i < querySet[0].size(); i++) {
        for (auto e:NodeSet[querySet[0][i]].categories)
        vector<PNode> pNodeList = FilterType(NodeSet, querySet[0][i], i,
                                             querySet, r);

        vector<PNode> pNodeListSpatial(pNodeList.begin(), pNodeList.end());
        if (!BASELINE) {
            sort(pNodeListSpatial.begin(), pNodeListSpatial.end(), comp2);
            sort(pNodeList.begin(), pNodeList.end(), comp);
        }
        oriList.push_back(pNodeList);
        oriListSpatial.push_back(pNodeListSpatial);
        possible *= pNodeList.size();
    }

    total_possible += possible;


    for (int i = 0; i < 10; i++)
        numMapping[i] = i;
    int ListSizes[10];
    for (int i = 0; i < oriList.size(); i++) {
        ListSizes[i] = oriList[i].size();
    }



    int oriListSize = oriList.size();

    double wz = 0;


    minLat.clear();
    minLon.clear();
    maxLat.clear();
    maxLon.clear();
    double totminLat = 1e18;
    double totminLon = 1e18;
    double totmaxLat = -1e18;
    double totmaxLon = -1e18;

    for (int i = 0;i < oriListSize; i++) {
        int oriListItemSize = oriList[i].size();
        minLat.push_back(1e18);
        minLon.push_back(1e18);
        maxLat.push_back(-1e18);
        maxLon.push_back(-1e18);
        wz += 1.0 / oriListItemSize;

        for (int j = 0; j < oriList[numMapping[i]].size(); j++) {
            minLat[i] = min(minLat[i], NodeSet[oriList[numMapping[i]][j].id].latitude);
            minLon[i] = min(minLon[i], NodeSet[oriList[numMapping[i]][j].id].longitude);
            maxLat[i] = max(maxLat[i], NodeSet[oriList[numMapping[i]][j].id].latitude);
            maxLon[i] = max(maxLon[i], NodeSet[oriList[numMapping[i]][j].id].longitude);
        }

        totminLat = min(totminLat, minLat[i]);
        totminLon = min(totminLon, minLon[i]);
        totmaxLat = max(totmaxLat, maxLat[i]);
        totmaxLon = max(totmaxLon, maxLon[i]);

    }

    double minmaxLat;
    double minmaxLon;

    for (int i = 0;i < oriListSize; i++) {

        minmaxLat = Spherical_distance(minLat[i],minLon[i],maxLat[i],minLon[i]);
        minmaxLon = Spherical_distance(minLat[i],minLon[i],minLat[i],maxLon[i]);


    }

    minmaxLat = Spherical_distance(totminLat,totminLon,totmaxLat,totminLon);
    minmaxLon = Spherical_distance(totminLat,totminLon,totminLat,totmaxLon);

    LatSplitNum = int(minmaxLat/100);
    LonSplitNum = int(minmaxLon/100);
    latLength = (totmaxLat - totminLat) / LatSplitNum;
    lonLength = (totmaxLon - totminLon) / LonSplitNum;





    extern double CostInPreSort;

    long start = clock();



    for (int i = 0; i < querySet.size(); i++) {
        vector<int> temp;
        temp.resize(oriListSize);
        for (int j = 0; j < oriListSize; j++) {
            temp[j] = querySet[i][numMapping[j]];
            //  temp[j] = querySet[i][j];
        }
        changedQuerySet.push_back(temp);
    }

    SpatialVector.clear();
    for (int i = 0; i < changedQuerySet.size(); i++) {
        int mag_sum = 0;
        for (int j = 0; j < changedQuerySet[i].size(); j++)
            mag_sum += changedQuerySet[i][j];
        mag_sum /= changedQuerySet[0].size();
        magnitude = mag_sum;
        SpatialVector.push_back(getSpatialVector(changedQuerySet[i]));
    }
    COUNT = 0;
    vector<int> combination(oriListSize);
    while (!que.empty()) {
        que.pop();
    }
    vector<double>  spatial_query_vec = getSpatialVector(changedQuerySet[0]);
    double length_limit = 0;
    for (auto wi : spatial_query_vec) {
        length_limit += wi * wi;
    }
    length_limit = sqrt(length_limit);
    splitDFS(NodeSet, oriList, changedQuerySet, length_limit, D, 0,  combination, que, 0, 1, totminLat, totmaxLat+1e-8, totminLon, totmaxLon+1e-8);


    int kcnt = 2;

    double UPbound = 1.0;




    vector<QueueNode> ans;
    while (!que.empty()) {
        ans.push_back(que.top());
        que.pop();
    }
    reverse(ans.begin(), ans.end());

    percentage += (possible - COUNT) * 1.0 / possible;

    return ans;
}

bool checkDist(double R, int a, int b) {
    return Spherical_distance(NodeSet[a].latitude, NodeSet[a].longitude,
                              NodeSet[b].latitude, NodeSet[a].longitude) <= R;
}

bool checkType(int a, int b) // b is the corresponding query node
{
    set<string>::iterator it;
    for (it = NodeSet[a].categories.begin(); it != NodeSet[a].categories.end();
         it++) {
        if (NodeSet[b].categories.find(*it) != NodeSet[b].categories.end())
            return true;
    }
    return false;
}


vector<vector<vector<int> > > ReadQuery(string path) {
    ifstream ifile;
    ifile.open(path.c_str());
    if (ifile.fail()) {
        cout << "file not found in ReadQuery" << endl;
    }
    int N;
    ifile >> N;
    vector<vector<vector<int> > > QuerySetList;
    RadioSet.clear();
    int temp;
    for (int i = 0; i < N; i++) {
        double R;
        int exampleNum;
        int exampleSize;
        ifile >> R >> exampleNum >> exampleSize;
        RadioSet.push_back(R);
        vector<vector<int> > QuerySet;
        for (int j = 0; j < exampleNum; j++) {
            vector<int> example;

            for (int k = 0; k < exampleSize; k++) {
                ifile >> temp;
                example.push_back(temp);
            }
            if (ONLY_ONE_EX == 1 && j == 0 || ONLY_ONE_EX != 1)
                QuerySet.push_back(example);
        }
        QuerySetList.push_back(QuerySet);
    }
    ifile.close();
    return QuerySetList;
}


int main(int argc, char* argv[]) {
//    std::srand ( unsigned ( std::time(0) ) );
    extern vector<Node> NodeSet;
    extern int K;

    string data_num = argv[1];
    string split_k = argv[2];
    int D = atoi(argv[3]);
    K = atoi(argv[4]);
    int case1 = atoi(argv[5]);
    int case2 = atoi(argv[6]);
    string query_type = argv[7];
    double alpha_v = atof(argv[8]);
    rate_p = atof(argv[9]);
    Error_bound_1 = atof(argv[10]);

    Read(NodeSet, "./gaode_poi_2018_value_"+data_num+".csv");
    for (auto e:NodeSet[0].categories)
        cout<<"e: "<<e<<endl;
    BuildIndex();

    string QueryFilePath = "gaode_query_";



    // a map from category to a list of node ids that contain this category
    int Q = 2;

    string path;

    if (query_type == "k_avg") {
        path = "./" + QueryFilePath + "num_" + data_num + "_type_num_3_split_k_avg.csv";
    }
    vector<vector<vector<int> > > QuerySet;

    QuerySet = ReadQuery(path);

//	case_study();
    int kSet[5] = { 1, 5, 10, 20, 50 };



    clock_t start;

    alpha = alpha_v;
    start = clock();
    CostInSpatialUb = 0;
    CostInPreSort = 0;
    CostInPreNear = 0;
    CostInSimMap = 0;
    percentage = 0;

    dfs_time = 0;
    int CASE = 100;
    int cas;

    for (cas = case1; cas < case2; cas++) {
        clock_t start1 = clock();


        vector<QueueNode> ans = getTopK(NodeSet, QuerySet[cas],
                                        RadioSet[cas], D, K);

        for (int i = 0; i < ans.size(); i++) {

            for (int j = 0; j < ans[i].id.size(); j++) {
                cout << ans[i].id[j] << "(" << NodeSet[ans[i].id[j]].latitude << "," << NodeSet[ans[i].id[j]].longitude << ")\t";
            }
            cout<<endl;

        }

        clock_t end1 = clock();


    }





    return 0;
}

