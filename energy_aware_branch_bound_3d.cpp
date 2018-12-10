#include<iostream>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include<vector>
#include<limits.h>
#include<string.h>
#include<cmath>
#include<sys/types.h>
#include<fcntl.h>
#include<sys/stat.h>
#include<cstdio>
#include<ctime>
using namespace std;
#define MAX_LINE 1024
#define link_length 200
char folder[100];
int dimension[3];
int ver;
char load_condition;
string folder_name;
string path;
string path1;
string path2;
string path3;
string rtable;
string topology;
string traff;
string rname;
string rem="rm -f ";
char apcgfilename[50];
bool parse_apcg(char * filename);
int ** appMatrix=NULL;
float ** icnEnergyMatrix=NULL;
int * app_rank_array=NULL;
int ** outgoing_bw_requirement=NULL;
int ** incoming_bw_requirement=NULL;
int * link_bw=NULL;
void build_icnEnergyMatrix();
int gTileNum,gLinkNum,gProcNum;
int g_edge_size;
int numrows;
int slice_num;
int s_size;
int a;
int snum;
int min_hit_threshold=40;
int pq_size=200;
int *** link_usage_matrix=NULL;  
vector<int> ** link_usage_list=NULL;
const float MAX_VALUE = INT_MAX - 100;
float MAX_PER_TRAN_COST = MAX_VALUE;
void initialize();
bool exist_locked_pe();
void bnb_map();
void randommapping();
void randommapping1();
void BBMClear();
void RANDOMMClear();
int counting=0;
void sorting_process();
void build_appMatrix();
void build_link_usage_matrix();
void delete_link_usage_matrix();
void delete_link_usage_list();
void generate_topology_IR_config();
void generate_link_length();
void file_exist(string rtable_filename1);
void create_files();
void remove_unwanted_files();
char type;
class traffic
{
public:
vector<int> load;
vector<int> dst;
vector<int> toVolume;
}; 
traffic t[1000];
typedef struct Position
{
    int row;
    int col;
    int slice_id;
} *pPosition, Position;
bool operator==(const Position & pos1, const Position &pos2);
//////////////////////////////////////////////////////////////////// Class Link//////////////////////////////////////////////////////////////////////////////////
typedef class Link{
  // friend ostream & operator<<(ostream &os, const Link & link);
    static int cnt;
    int id;
    Position fromTile;
    Position toTile;
    int fromTileId;
    int toTileId;
    static int count;
    static int slicenum;
    static int lnk;
    static int lnk1;
    static int k;
    static bool status;
    static int record;
    static float powerCoeff;
    public:
    Link();
    float Cost() {return powerCoeff;}
    const Position & FromTilePos();
    const Position & ToTilePos();
    int FromTile() { return fromTileId;}
    int ToTile() {return toTileId;}
    int const GetCnt();
    int GetId() {return id;}
    }*pLink, Link;
///////////////////////////////////////////////////////////Class Link End//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////Class Tile//////////////////////////////////////////////////////////////////////////////////////////////////
typedef class Tile
{
public:
	static int cnt;
	int id;
	static int count;
	static int slicenum;
	Position pos;
	int goLinkNum;
	int comeLinkNum;
	vector<pLink> goLink;
	vector<pLink> comeLink;
	float powerCoeff;
	int procId;
	class Router{
        Tile * host_tile;               
        int **routing_table;      //routing_table[i][j] shows how to which link to send the 
        //from tile i to tile j. If it's -2, means not reachable.
        //if it's -1, means the destination is the current tile.
        bool generate_xyz_routing_table();
    public:
        Router();
        bool initialize(Tile * host);
        int route_to_link(int src_tile, int dst_tile) const;
        int set_routing_entry(int src_tile, int dst_tile, int link_id);
        ~Router();
        } router;
	float Cost()
	     {
	         return powerCoeff;
	     }
	Tile();
	int GetId() const;
	int GetGoLinkNum() const;
	int GetComeLinkNum() const;
	Position GetPosition() const;
	pLink GoLink(int i);
	pLink ComeLink(int i);
	void AttachLink(pLink gLink);
        bool initialize_router();
        int RouteToLink(int srcId, int dstId) const;
        int RouteToLink(const Tile& srcTile, const Tile& dstTile) const;
        int RouteToLink(int srcId, int dstId, int linkId);
        int set_routing_entry(int src_tile, int dst_tile, int link_id);
	//friend ostream & operator<<(ostream & os, const Tile & tile);
	friend bool operator==(const Position & pos1, const Position &pos2);
	~Tile();
}*pTile, Tile;
//////////////////////////////////////////////Class Tile End///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////Class Router/////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef class Router{
    Tile * host_tile;               
    int **routingTable;      //routingTable[i][j] shows how to which link to send the 
    //from tile i to tile j. If it's -2, means not reachable.
    //if it's -1, means the destination is the current tile.

    bool generate_xyz_routing_table();
   
public:
    Router();
    bool initialize(pTile host);
    int route_to_link(int src_tile, int dst_tile) const;
    int set_routing_entry(int src_tile, int dst_tile, int link_id);
    ~Router();
} *pRouter, Router;
//////////////////////////////////////////////Class Router End////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////Class Mapping Node////////////////////////////////////////////////////////////////////////////////////////////////////
class MappingNode
{
public:
static int cnt;
void initialize();
friend void bnb_map();
friend void randommapping();
friend void randommapping1();
     bool illegal;    // It is an illegal node if it violates the spec
    // constructor will init this
     int stage;       // How many processes have been mapped
     int *mappingSequency;
     bool *tileOccupancyTable;
     long cost;
     long llc; //lower limit cost
     long ulc; //upper limit cost
     bool occupancyTableReady;
     void GreedyMapping();
     long llc();
     long ulc();
     long lowestUnitCost(int tileId);
     long lowestUnmappedUnitCost();
     void MapNode(int procId, float goodRow, float goodCol, float goodSlice_id);
     int *link_BW_usage;
     bool illegal_child_mapping;
     class MappingNode* next;
     friend class PQueue;
 bool Expandable(int);		
bool fixed_verify_BW_usage();
     int bestulcCandidate();
     int bestCostCandidate();
MappingNode()
{
} 
bool is_illegal(void) const
 {
        return illegal;
 }
int mapToTile(int i)
    {
         return mappingSequency[i];
    }
     MappingNode(int tileId);
     MappingNode(const MappingNode& parent, int tileId, bool calcBound=true);
	
};
int MappingNode::cnt;
typedef class MappingNode MappingNode, *pMappingNode;
//////////////////////////////////////////////////////////Class Mapping Node End///////////////////////////////////////////////////////////////////////
void BBMMap(pMappingNode bestMapping);
void randommap(pMappingNode randomMapping);
void generate_rtable(pMappingNode bestMapping);
////////////////////////////////////////////////////////////Class Priority Queue/////////////////////////////////////////////////////////////////////////
 // A class for priority queue.
 // TODO: use heapsort data structure for speeding up
 typedef class PQueue
 {
 private:
     int length;
     pMappingNode head;
     int minCost;
     int minllc;
     int minulc;
 public:
     PQueue();
     ~PQueue();
     int Length()
     {
         return length;
     }
     bool empty();
     void insert(pMappingNode node);
     pMappingNode next();
pMappingNode traverse();
 } PQueue, *pPQueue;
///////////////////////////////////////////////////////////Class Priority Queue End//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////Class Process//////////////////////////////////////////////////////////////////////////////////////////
class Process {
    static int proc_num;              // the number of processes in the application
 public:
    int id;                    // the unique id of that process, from 0 -> app_size-1
    int *toComm; 
    int *fromComm;
    int *outgoing_bw_requirement;    // the bandwidth requirement of the out-going traffic
    int *incoming_bw_requirement;
    int rank;
    int tileId;
    int totalCommVol;
    int mapToTile(int tid)
 	{
         tileId = tid;
         return tid;
 	}
    Process();	
friend void initialize();
 };
int Process::proc_num; 
typedef class Process *pProcess, Process;
/////////////////////////////class Process End//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<pProcess> gProcess;
Process::Process() {
    id = proc_num++;
    toComm = NULL;
    fromComm = NULL;
    outgoing_bw_requirement = NULL;
    incoming_bw_requirement = NULL;
    toComm = new int[gProcNum];
    fromComm = new int[gProcNum];
    outgoing_bw_requirement = new int[gProcNum];
    incoming_bw_requirement = new int[gProcNum];
    for (int i=0; i<gProcNum; i++) 
{
        toComm[i] = 0;
        fromComm[i] = 0;
        outgoing_bw_requirement[i] = NULL;
        incoming_bw_requirement[i] = NULL;
}
}
float Link :: powerCoeff=0.00007;
pLink gLink = NULL;
pTile gTile = NULL;
int Tile::cnt=0;
int Tile::count;
int Tile::slicenum;
Tile::Tile()
{
if(slice_num==1)
powerCoeff=0.52;
else if(slice_num>=2)
powerCoeff=0.54;
id=cnt++;
if(count==s_size)
{
count=0;
slicenum++;
}
//cout<<"current id = "<<id<<endl;
pos.row = count / g_edge_size;
pos.col = count % g_edge_size;
pos.slice_id=slicenum;
//cout<<"tile-"<<id<<" = "<<pos.row<<pos.col<<pos.slice_id<<endl;
count++;
}
int Link::cnt;
int Link::count;
int Link::slicenum;
int Link::lnk;
int Link::lnk1;
int Link::k;
int Link::record;
bool Link::status=false;
Link::Link() { 
id = cnt++;
int x=((2 * (g_edge_size-1) * numrows) +((numrows-1) * 2 * g_edge_size));
if(count==x && record < (2*s_size) && slicenum<slice_num)
{
//cout<<"record = "<<record<<endl;
//cout<<"lnk = "<<lnk<<endl;
int a=abs((s_size*(slicenum+1))-s_size);
//cout<<"a = "<<a<<endl;
if(k<g_edge_size && status==false)
{
int localid=a+lnk;
int src=localid;
int dst=localid+s_size;
//cout<<"src="<<src<<endl;
//cout<<"dst="<<dst<<endl;
fromTile.row=(src%s_size)/g_edge_size;
toTile.row=fromTile.row;
fromTile.col=(src%s_size)%g_edge_size;
toTile.col=(dst%s_size)%g_edge_size;
fromTile.slice_id=slicenum;
toTile.slice_id=slicenum+1;
k++;
//cout<<"K = "<<k<<" ";
if(k==g_edge_size)
status=true;
//cout<<"and status = "<<status<<endl;
lnk++;
}
else if(k<=g_edge_size && status==true)
{
//cout<<"hello"<<endl;
int localid=a+lnk1;
//int localid=(localid)%(g_edge_size-1);
//cout<<"localid = "<<localid<<endl;
int src=localid+s_size;
int dst=localid;
//cout<<"src="<<src<<endl;
//cout<<"dst="<<dst<<endl;
toTile.row=(dst%s_size)/g_edge_size;
fromTile.row=toTile.row;
fromTile.col=(src%s_size)%g_edge_size;
toTile.col=(dst%s_size)%g_edge_size;
fromTile.slice_id=slicenum+1;
toTile.slice_id=slicenum;
k--;
//cout<<"k = "<<k<<endl;
if(k==0)
status=false;
lnk1++;
}
record++;
if(record==2*s_size)
{
count=0;
slicenum++;
record=0;
lnk=0;
lnk1=0;
}
}
else
{
//cout<<"Count = "<<count<<endl; 
//cout<<"Slicenum = "<<slicenum<<endl; 
    //There are totally 2*(g_edge_size-1)*g_edge_size*2 links. The first half links are horizontal
    //the second half links are veritcal links. 
    if (count < ((2 * (g_edge_size-1) * numrows))) {
        fromTile.row = count/(2*(g_edge_size-1));
        toTile.row = count/(2*(g_edge_size-1));
        fromTile.slice_id=slicenum;
	toTile.slice_id=slicenum;
        int localId = count%(2*(g_edge_size-1));
        if (localId < (g_edge_size-1)) {
            //from west to east
            fromTile.col = localId;
            toTile.col = localId + 1;	
        }
        else {
            //from east to west
            localId = localId - (g_edge_size-1);
            fromTile.col = localId + 1;
            toTile.col = localId;
        }
    }
    else {
        int localId = count-((2 * (g_edge_size-1) * numrows));
        fromTile.col = localId/(2*(numrows-1));
        toTile.col = localId/(2*(numrows-1));
	fromTile.slice_id=slicenum;
	toTile.slice_id=slicenum;
        localId = localId%(2*(numrows-1));
        if (localId < (numrows-1)) {
            //from south to north
            fromTile.row = localId;
            toTile.row = localId + 1;
        }
        else {
            //from north to south
            localId = localId - (numrows-1);
            fromTile.row = localId + 1;
            toTile.row = localId;
        }
    }
count++;   
}
//For mesh, all the links have the same power coefficency
//cout<<"Link "<< id <<" created between tile-";
//cout<<fromTile.row<<fromTile.col<<fromTile.slice_id;
//cout<<" and tile-"<<toTile.row<<toTile.col<<toTile.slice_id<<endl;
fromTileId = fromTile.row*g_edge_size + fromTile.col+(s_size * fromTile.slice_id);
toTileId = toTile.row*g_edge_size + toTile.col+(s_size * toTile.slice_id);
//cout<<"Link "<< id <<" created between tile-"<< fromTileId<<" and tile-"<<toTileId<<endl;
}
int Tile::GetId() const {
    return id;
}

int Tile::GetGoLinkNum() const {
    return goLinkNum;
}

int Tile::GetComeLinkNum() const {
    return comeLinkNum;
}

pLink Tile::GoLink(int i) {
    return goLink[i];
}

pLink Tile::ComeLink(int i) {
    return comeLink[i];
}
Tile::~Tile() {
    goLink.clear();
    comeLink.clear();
}
Tile::Router::Router() {
    routing_table = new int*[gTileNum];
    for(int i=0; i<gTileNum; i++) 
        routing_table[i] = new int[gTileNum];
    for(int i=0; i<gTileNum; i++) 
        for(int j=0; j<gTileNum; j++) 
            routing_table[i][j] = -2;
}

Tile::Router::~Router() {
    if(routing_table) {
        for(int i=0; i<gTileNum; i++) 
            delete []routing_table[i];
        delete []routing_table;
    }
}
void Tile::AttachLink(pLink gLink) {
    goLinkNum = comeLinkNum = 0;
    for (int i=0; i<gLink[0].GetCnt(); i++) {
//cout<<"gLink[0].GetCnt()"<<gLink[0].GetCnt()<<endl;
        if (gLink[i].FromTilePos() == pos) {
            goLink.push_back(&gLink[i]);
 cout<< gLink[i].GetId()<< " ";
            goLinkNum++;
        }
        if (gLink[i].ToTilePos() == pos) {
           cout<< gLink[i].GetId()<< " ";
            comeLink.push_back(&gLink[i]);
            comeLinkNum++;
        }
    }
}
Position Tile::GetPosition() const
{
    return pos;
}
int const Link::GetCnt() {
    return cnt;
}

const Position & Link::FromTilePos() {
    return fromTile;
}

const Position & Link::ToTilePos() {
    return toTile;
}
bool Tile::initialize_router()
{
    return router.initialize(this);
}

int Tile::RouteToLink(int srcId, int dstId) const {
    return router.route_to_link(srcId, dstId);
}

int Tile::RouteToLink(const Tile& srcTile, const Tile& dstTile) const {
    return router.route_to_link(srcTile.id, dstTile.id);
}

//This method is used to program the routing table
int Tile::RouteToLink(int srcId, int dstId, int linkId) {
    return router.set_routing_entry(srcId, dstId, linkId);
}
bool operator==(const Position & pos1, const Position & pos2) 
{
    return (pos1.row==pos2.row && pos1.col==pos2.col && pos1.slice_id==pos2.slice_id);
}
bool Tile::Router::initialize(pTile host) 
{
    host_tile = host;
    return generate_xyz_routing_table();
}

///////////////////////////////////////XYZ routing/////////////////////////////////////////////////////

bool Tile::Router::generate_xyz_routing_table() {
cout<<"Generating routing table"<<endl;
    for(int i=0; i<gTileNum; i++) 
        for(int j=0; j<gTileNum; j++) 
            routing_table[i][j] = -2;
  
    for(int dstTile=0; dstTile<gTileNum; dstTile++) {
        if(dstTile == host_tile->id) {                     //deliver to me
            routing_table[0][dstTile] = -1;
            continue;
        }

        //check out the dst Tile's position first
        Position dstPos;
        dstPos.row = (dstTile%s_size)/g_edge_size;
        dstPos.col = (dstTile%s_size)%g_edge_size;
        dstPos.slice_id=dstTile/s_size;
        Position nextStep = host_tile->pos;
//cout<<"Source tile is tile-"<<host_tile->pos.row<<host_tile->pos.col<<host_tile->pos.slice_id<<" and "<<"destination tile is tile-"<<dstPos.row<<dstPos.col<<dstPos.slice_id<<endl;
        if(dstPos.col != host_tile->pos.col) {            //We should go horizontally
            if(host_tile->pos.col > dstPos.col) 
                nextStep.col--;
            else 
                nextStep.col++;
//cout<<"Next tile is tile-"<<nextStep.row<<nextStep.col<<nextStep.slice_id<<endl;
        }
        else if(dstPos.col == host_tile->pos.col && dstPos.row != host_tile->pos.row)
{            //We should go vertically
            if(host_tile->pos.row > dstPos.row) 
                nextStep.row--;
            else 
                nextStep.row++;
//cout<<"Next tile is tile-"<<nextStep.row<<nextStep.col<<nextStep.slice_id<<endl;
        }
else if(dstPos.row==host_tile->pos.row && dstPos.col==host_tile->pos.col && dstPos.slice_id!=host_tile->pos.slice_id)
{
if(host_tile->pos.slice_id > dstPos.slice_id) 
                nextStep.slice_id--;
            else 
                nextStep.slice_id++;
//cout<<"Next tile is tile-"<<nextStep.row<<nextStep.col<<nextStep.slice_id<<endl;
}

        int i=0;
        for(i=0; i<host_tile->goLinkNum; i++) {
            pLink pL = host_tile->goLink[i];
            if(pL->ToTilePos() == nextStep) {
                routing_table[0][dstTile] = pL->GetId();
                break;
            }
        }
    }

    //Duplicate this routing row to the other routing rows.
    for(int i=1; i<gTileNum; i++) 
        for(int j=0; j<gTileNum; j++) 
            routing_table[i][j] = routing_table[0][j];

    return true;
}
int Tile::Router::route_to_link(int src_tile, int dst_tile) const {
    return routing_table[src_tile][dst_tile];
}

int Tile::Router::set_routing_entry(int src_tile, int dst_tile, int link_id) {
    routing_table[src_tile][dst_tile] = link_id;
    return link_id;
}
MappingNode::MappingNode(const MappingNode& parent, int tileId, bool calcBound)
 {
     illegal = false;
     tileOccupancyTable = NULL;
     mappingSequency = NULL;
     link_BW_usage = NULL;
     occupancyTableReady = false;
     llc = -1;
     cnt ++;
     mappingSequency = new int[gProcNum];
     for (int i=0; i<gProcNum; i++)
         mappingSequency[i] = -1;

     stage = parent.stage;
 

     int proc1 = app_rank_array[stage];
     /*if (gProcess[proc1]->is_locked() && gProcess[proc1]->lock_to != tileId)
     {
         illegal = true;
         return;
     }*/
 
     llc = parent.llc;
     ulc = parent.ulc;
 
     // Copy the parent's partial mapping
     memcpy(mappingSequency, parent.mappingSequency, gProcNum*sizeof(int));
 
         link_BW_usage = new int[gLinkNum];
         memcpy(link_BW_usage, parent.link_BW_usage, sizeof(int)*gLinkNum);
     
 
     //Map the next process to tile tileId
     mappingSequency[stage] = tileId;
     next = NULL;
     cost = parent.cost;

/*char str_id[5];
	ofstream pout;
	string traffic_filename = string("output1/node-");
	// create new traffic log file
	for(int i=0; i<=stage; i++) {
		sprintf(str_id, "%d", mappingSequency[i]);
		cout<<str_id;
		traffic_filename+=string(str_id);
               }
		pout.open(traffic_filename.c_str());
		pout.close(); */
     for (int i=0; i<stage; i++)
     {
         int tile1 = tileId;
         int tile2 = mappingSequency[i];
         float this_tran_cost = appMatrix[i][stage];
         this_tran_cost = this_tran_cost * icnEnergyMatrix[tile1][tile2];
         cost += this_tran_cost;
         if (this_tran_cost > MAX_PER_TRAN_COST)
         {
             illegal = true;
             return;
         }
     }
// cout<<"Cost between tileid-"<<tileId<<" and mappingSequecny node-";
//for(int i=0;i<stage;i++)
//cout<<mappingSequency[i];
//cout<<" = "<<cost<<endl; 
     
         for (int i=0; i<stage; i++)
         {
             int tile1 = tileId;
             int tile2 = mappingSequency[i];
             int proc1 = app_rank_array[stage];
             int proc2 = app_rank_array[i];
             if (outgoing_bw_requirement[proc1][proc2])
             {
                 for (unsigned int i=0; i<link_usage_list[tile1][tile2].size(); i++)
                 {
                     int link_id = link_usage_list[tile1][tile2][i];
                     link_BW_usage[link_id] += outgoing_bw_requirement[proc1][proc2];
                     if (link_BW_usage[link_id] > link_bw[link_id])
                     {
                         cost = MAX_VALUE+1;
                         illegal = true;
                         return;
                     }
                 }
             }
             if (incoming_bw_requirement[proc2][proc1])
             {
                 for (unsigned int i=0; i<link_usage_list[tile2][tile1].size(); i++)
                 {
                     int link_id = link_usage_list[tile2][tile1][i];
                     link_BW_usage[link_id] += incoming_bw_requirement[proc2][proc1];
                     if (link_BW_usage[link_id] > link_bw[link_id])
                     {
                         cost = MAX_VALUE+1;
                         illegal = true;
                         return;
                     }
                 }
             }
         }
    
 
     stage ++;
 
     if (!calcBound)
         return;
 
     tileOccupancyTable = new bool[gTileNum];
     for (int i=0; i<gTileNum; i++)
         tileOccupancyTable[i] = false;
 
     llc = llc();
     ulc = ulc();
//cout<<"Cost = "<<cost<<endl;
 }
 
 MappingNode::MappingNode(int tileId)
{

	// initialize random number generator
	/*char str_id[5];
	ofstream pout;
	
	// create new traffic log file
	//for(int i = 0; i < num_tiles; i++) {
		sprintf(str_id, "%d", tileId);
		cout<<str_id;
		string traffic_filename = string("output1/node-") + string(str_id);
		pout.open(traffic_filename.c_str());
		pout.close();*/
     illegal = false;
     tileOccupancyTable = NULL;
     mappingSequency = NULL;
     link_BW_usage = NULL;
     occupancyTableReady = false;
     llc = -1;
 
     cnt ++;
     mappingSequency = new int[gProcNum];
     for (int i=0; i<gProcNum; i++)
         mappingSequency[i] = -1;
 
     stage = 1;
     mappingSequency[0] = tileId;
     next = NULL;
     cost = 0;
 
     int proc1 = app_rank_array[0];
    /* if (gProcess[proc1]->is_locked() && gProcess[proc1]->lock_to != tileId)
     {
         illegal = true;
         return;
     }*/
 
 
     tileOccupancyTable = new bool[gTileNum];
     for (int i=0; i<gTileNum; i++)
         tileOccupancyTable[i] = false;
 
        link_BW_usage = new int[gLinkNum];
        for (int i=0; i<gLinkNum; i++)
             link_BW_usage[i] = 0;
 
     llc = llc();
     ulc = ulc();
//cout<<"Cost = "<<cost<<endl;
 }
void initialize()
{
//cout<<"creates and initializes icnEnergyMatrix and parse the given apcg input file"<<endl; 
app_rank_array=new int[gTileNum];
build_icnEnergyMatrix();
link_bw=new int[gLinkNum];
for(int i=0;i<gLinkNum;i++)
link_bw[i]=1000000;
outgoing_bw_requirement = new int*[gProcNum];
for (int i=0; i<gProcNum; i++)
outgoing_bw_requirement[i] = new int[gProcNum];
incoming_bw_requirement = new int*[gProcNum];
for (int i=0; i<gProcNum; i++)
incoming_bw_requirement[i] = new int[gProcNum];
 for (int i=0; i<gProcNum; i++)
    {
           pProcess p = new Process();
           gProcess.push_back(p);
    }            
//parse_apcg("apcg-8.txt");
//parse_apcg("apcg-36.txt");
parse_apcg(apcgfilename);
//parse_apcg("apcg-9.txt");
//parse_apcg("apcg-12.txt");
//parse_apcg("apcg-5.txt");
build_link_usage_matrix();
}
void build_link_usage_matrix() {
    if (link_usage_matrix) 
        delete_link_usage_matrix();
    // Allocate the space for the link usage table
    link_usage_matrix = new int**[gTileNum];
    for (int i=0; i<gTileNum; i++) {
        link_usage_matrix[i] = new int*[gTileNum];
        for (int j=0; j<gTileNum; j++) 
            link_usage_matrix[i][j] = new int[gLinkNum];
    }

    for (int i=0; i<gTileNum; i++) 
        for (int j=0; j<gTileNum; j++) 
            for (int k=0; k<gTileNum; k++) 
                link_usage_matrix[i][j][k] = 0;

    // Setting up the link usage matrix
    for (int src_id=0; src_id<gTileNum; src_id++) {
        for (int dst_id=0; dst_id<gTileNum; dst_id++) {
            if (src_id == dst_id)
                continue;
            pTile current_tile = &gTile[src_id];
            while (current_tile->GetId() != dst_id) {
                int link_id = current_tile->RouteToLink(src_id, dst_id);
                pLink pL = &gLink[link_id];
                link_usage_matrix[src_id][dst_id][link_id] = 1;
                current_tile = &gTile[pL->ToTile()];
            }
        }
    }
  
    if (link_usage_list) 
        delete_link_usage_list();

    //Now build the g_link_usage_list
    link_usage_list = new vector<int>*[gTileNum];
    for (int i=0; i<gTileNum; i++) 
        link_usage_list[i] = new vector<int>[gTileNum];
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {
            link_usage_list[src][dst].clear();
            if (src==dst) 
                continue;
            for (int link_id=0; link_id<gLinkNum; link_id++) {
                if (link_usage_matrix[src][dst][link_id]) 
                    link_usage_list[src][dst].push_back(link_id);
            }
        }
    }
}
long MappingNode::llc()
 {
//cout<<mappingSequency[0];
//ofstream pout;
//pout.open("outputnew.txt",ios::app);
//pout<<"stage="<<stage<<endl;
//for(int i=0;i<16;i++)
//{
//for(int j=0;j<16;j++)
//{
//pout<<appMatrix[i][j]<<endl;
//}
//}
char str_id[500];
	ofstream pout;
	string traffic_filename = string("output111/node-");
	// create new traffic log file
	for(int i=0; i<stage; i++) {
		sprintf(str_id, "%d", mappingSequency[i]);
		traffic_filename+=string(str_id);
               }
		pout.open(traffic_filename.c_str(), ios::app);
		
    for (int i=0; i<gTileNum; i++)         
	tileOccupancyTable[i] = false;
	for (int i=0; i<stage; i++)
tileOccupancyTable[mappingSequency[i]] = true;
    occupancyTableReady = true;
llc = cost;
   for (int i=0; i<stage; i++)
    {
         for (int j=stage; j<gProcNum; j++)
         {
//pout<<appMatrix[i][j]<<endl;     
        if (appMatrix[i][j]==0)
                 continue;
             else
                 llc += appMatrix[i][j]*lowestUnitCost(mappingSequency[i]);
	pout<<"llc for i="<<i<<" and j="<<j<<" is :"<<" "<<llc<<endl;
         }
     }
     //Now add the cost of the communication among all the un-mapped nodes
     int vol = 0;
     for (int i=stage; i<gProcNum; i++)
     {
         for (int j=i+1; j<gProcNum; j++)
	{
             vol += appMatrix[i][j];
		pout<<"Volume for i= "<<i<<" and j= "<<j <<" is : "<<vol<<endl;
}
     }
     llc += vol * lowestUnmappedUnitCost();
pout<<"Total llc ater adding the cost of the communication among all the un-mapped nodes is: "<<llc<<endl;
pout.close();
//cout<<endl<<"llc of node-";
//for(int i=0;i<stage;i++)
//cout<<mappingSequency[i];
//cout<<" = "<<llc<<endl;

     return llc;
 }

long MappingNode::lowestUnitCost(int tileId)
 {
     long min = 50000;

     for (int i=0; i<gTileNum; i++)
     {
         if (i==tileId)
             continue;
         if (tileOccupancyTable[i])
             continue;
         if (icnEnergyMatrix[tileId][i] < min)
             min = icnEnergyMatrix[tileId][i];
     }
//pout<<gTileNum<<endl;
//pout<<gProcNum<<endl;
     return min;
 } 
 //This function returns the lowest cost between anytwo unoccupied tiles
 long MappingNode:: lowestUnmappedUnitCost()
 {
     long min = 50000;
//pout<<gTileNum<<endl;
//pout<<gProcNum<<endl;
//for (int i=0; i<gTileNum; i++)
  //   {
    //     pout<<tileOccupancyTable[i]<<" "<<i<<endl;
//}
   for (int i=0; i<gTileNum; i++)
     {
         if (tileOccupancyTable[i])
             continue;
         for (int j=i+1; j<gTileNum; j++)
         {
             if (tileOccupancyTable[j])
                 continue;
             if (icnEnergyMatrix[i][j]<min)
                 min = icnEnergyMatrix[i][j];
         }
     }
     return min;
} 
/////////////////////Upper Bound////////////////////////////////////////////
//This calculate the upper bound cost of the this partial mapping
 //in the current mapping
long MappingNode::ulc()
 {
//ofstream pout;
//pout.open("outputnew.txt",ios::app);
char str_id[500];
	ofstream pout;
	string traffic_filename = string("output111/node-");


	// create new traffic log file
	for(int i=0; i<stage; i++) {
		sprintf(str_id, "%d", mappingSequency[i]);
		
		traffic_filename+=string(str_id);
               }
		pout.open(traffic_filename.c_str(), ios::app);
		
pout<<"INSIDE ulc"<< endl<<"occupancyTableReady Status: "<< occupancyTableReady<<endl;
//cout<<"ulc of node having mappingSequency node-";
for (int i=0; i<stage; i++)
{
//cout<<mappingSequency[i];
 pout<<"tileOccupancyTable["<<mappingSequency[i]<<"] status = "<<tileOccupancyTable[mappingSequency[i]]<<endl;
}
cout<<endl;
     if (!occupancyTableReady)
     {
         for (int i=0; i<gTileNum; i++)
             tileOccupancyTable[i] = false;
         for (int i=0; i<stage; i++)
             tileOccupancyTable[mappingSequency[i]] = true;
     }

     GreedyMapping();
for(int i=0;i<gProcNum;i++)
pout<<"mappingSequency["<<i<<"] = "<< mappingSequency[i]<<endl;
     ulc = cost;
illegal_child_mapping=false;
 if (!fixed_verify_BW_usage())
     {
       illegal_child_mapping = true;
       ulc = MAX_VALUE;
        return ulc;
     }
//cout<<"ulc = "<<ulc<<endl;
 for (int i=0; i<stage; i++)
  {
         int tile1 = mappingSequency[i];
         for (int j=stage; j<gProcNum; j++)
         {
             int tile2 = mappingSequency[j];
pout<<"for tile-"<<i<<" and tile-"<<j<<"ulc+=appMatrix["<<i<<"]["<<j<<"] * icnEnergyMatrix["<<tile1<<"]["<<tile2<<"]=";
             ulc += appMatrix[i][j] * icnEnergyMatrix[tile1][tile2];
pout<<" "<<ulc<<endl;
         }
     }
pout<<"Total ulc1 = "<<ulc<<endl;
     for (int i=stage; i<gProcNum; i++)
     {
         int tile1 = mappingSequency[i];
         for (int j=i+1; j<gProcNum; j++)
         {
             int tile2 = mappingSequency[j];
pout<<"for tile-"<<i<<" and tile-"<<j<<"ulc+=appMatrix["<<i<<"]["<<j<<"] * icnEnergyMatrix["<<tile1<<"]["<<tile2<<"]=";
             ulc += appMatrix[i][j] * icnEnergyMatrix[tile1][tile2];
pout<<" "<<ulc<<endl;
         }
     }
pout<<"Total ulc2 = "<<ulc<<endl;
pout<<"Upper Bound is:"<<ulc<<endl;
 pout.close();
//cout<<"ulc of node-";
//for(int i=0;i<stage;i++)
//cout<<mappingSequency[i];
//cout<<" = "<<ulc<<endl;

     return ulc;
}
//Map the other unmapped process node using greedy mapping
 void MappingNode::GreedyMapping()
 {
char str_id[500];
	ofstream pout;
	string traffic_filename = string("output111/node-");


	// create new traffic log file
	for(int i=0; i<stage; i++) {
		sprintf(str_id, "%d", mappingSequency[i]);
		
		traffic_filename+=string(str_id);
               }
		pout.open(traffic_filename.c_str(), ios::app);
pout<<"INSIDE GREEDY MAPPING"<<endl;
     for (int i=stage; i<gProcNum; i++)
     {
	 int cnt=0;
         int sumRow = 0;
         int sumCol = 0;
         int sumSlice_id = 0;
         int vol = 0;
         for (int j=0; j<i; j++)
         { pout<<"for i = "<<i<<" and j = "<<j<<endl<<"appMatrix["<<i<<"]["<<j<<"] = "<<appMatrix[i][j]<<endl;
             if (appMatrix[i][j]==0)
                continue;
             int tileId = mappingSequency[j];
             int row = (tileId%s_size)/g_edge_size;
             int col = (tileId%s_size)%g_edge_size;
	     int Slice_id=tileId/s_size;
             sumRow += appMatrix[i][j]*row;
             sumCol += appMatrix[i][j]*col;
	     sumSlice_id += appMatrix[i][j]*Slice_id;
             vol += appMatrix[i][j];
	     cnt++;
pout<<"row = "<< row<<" col= "<< col<<" Slice_id= "<< Slice_id<<endl;
pout<<"sumRow = "<< sumRow<<" sumCol= "<< sumCol<<" sumSlice_id= "<< sumSlice_id<<endl;
pout<<"volume for i "<<i<<" and j= "<<j<<" is : "<< vol<<endl;
         }
         //This is somehow the ideal position
         float myRow, myCol,mySlice_id;
         if (vol==0)
         {
             myRow = -1;
             myCol = -1;
	     mySlice_id = -1;		

         }
         else
         {
             myRow = ((float) sumRow)/vol;
             myCol = ((float) sumCol)/vol;
             mySlice_id = ((float) sumSlice_id)/vol;
         }
/*int vol1=vol/cnt;
     int a=0;
	 for (int j=0; j<i; j++)
		{
			if (appMatrix[i][j]==0)
                	continue;
			if(vol1==appMatrix[i][j])
			a++;			
		}
if(cnt==a && a>1)
{
pout<<"Calling mapnode(i/procId, myRow, myCol) as : mapnode1( "<<i<<" ,"<<myRow<<" ,"<<myCol<<")"<<endl; 
MapNode1(i,myRow,myCol);
}
else
{*/
pout<<"Calling mapnode(i/procId, myRow, myCol, mySlice_id) as : mapnode( "<<i<<" ,"<<myRow<<" ,"<<myCol<<" ,"<<mySlice_id<<")"<<endl; 
         MapNode(i, myRow, myCol, mySlice_id);
//}
}
pout.close();
 }
 
 //Try to map the node to an unoccupied tile which is cloest to
 //the tile(goodRow, goodCol)
 void MappingNode::MapNode(int procId, float goodRow, float goodCol, float goodSlice_id)
 {
//ofstream pout;
//pout.open("outputnew.txt",ios::app);
char str_id[500];
	ofstream pout;
	string traffic_filename = string("output111/node-");
	// create new traffic log file

	for(int i=0; i<stage; i++) {
		sprintf(str_id, "%d", mappingSequency[i]);
		
		traffic_filename+=string(str_id);
               }
		pout.open(traffic_filename.c_str(), ios::app);
		 
     float minDist = 10000;
     int bestId = -1;
     for (int i=0; i<gTileNum; i++)
     {
         if (tileOccupancyTable[i])
{
      
pout<<"TileId "<< i<<" is occupied"<<endl;
       continue;
}
         if (goodRow<0)
         {
             bestId = i;
             break;
         }
         int row = (i%s_size)/g_edge_size;
         int col = (i%s_size)%g_edge_size;
	 int Slice_id=i/s_size;
pout<<"Inside mapnode function for tileId: "<<i<<" -> row = "<<row<<" col = "<<col<<" slice_id = "<<Slice_id<<endl;
pout<<"fabs(goodRow-row) ="<<fabs(goodRow-row)<<", fabs(goodCol-col) = "<<fabs(goodCol-col)<<"and fabs(goodSlice_id-Slice_id) = "<<fabs(goodSlice_id-Slice_id)<<endl;
         float dist = fabs(goodRow-row) + fabs(goodCol-col) + fabs(goodSlice_id-Slice_id);
pout<<"dist =  fabs(goodRow-row) + fabs(goodCol-col) : "<<dist<<endl;
         if (dist<minDist)
         {
             minDist = dist;
             bestId = i;
         }
pout<<"Current minimum distance is: "<<minDist<<endl;
pout<<"Current bestId is: "<<bestId<<endl;
     }
     mappingSequency[procId] = bestId;
     tileOccupancyTable[bestId] = true;
pout<<"mappingSequency["<<procId<<"] = "<<bestId<<endl;
pout<<"tileOccupancyTable["<<bestId<<"] = true"<<endl;
pout.close();
}
void delete_link_usage_list()
{
     if (link_usage_list)
     {
         for (int i=0; i<gTileNum; i++)
             delete []link_usage_list[i];
         delete []link_usage_list;
    }
 }
 
void delete_link_usage_matrix()
 {
     //free g_link_usage_matrix
     if (link_usage_matrix)
     {
         for (int i=0; i<gTileNum; i++)
             for (int j=0; j<gTileNum; j++)
                 delete []link_usage_matrix[i][j];
 
         for (int i=0; i<gTileNum; i++)
             delete []link_usage_matrix[i];
         delete []link_usage_matrix;
     }
 }
 PQueue::PQueue()
 {
     length = 0;
     head = NULL;
 }
 
 PQueue::~PQueue()
 {
 }
 
 bool PQueue::empty()
 {
     if (length==0)
         return true;
     else
         return false;
 }
 
 void PQueue::insert(pMappingNode node)
 {
     // here we should insert the node at the position which
     // is decided by the cost of the node
     if (length == 0)
     {
         head = node;
         node->next = NULL;
         length++;
         return;
     }
     pMappingNode parentNode = NULL;
     pMappingNode curNode = head;
     int i=0;
     for (i=0; i<length; i++)
     {
         if (curNode->cost > node->cost)
             //if (curNode->ulc > node->ulc)
             break;
         parentNode = curNode;
         curNode = curNode->next;
     }
     if (parentNode == NULL)
     {
         pMappingNode oldHead = head;
         head = node;
         node->next = oldHead;
         length++;
         return;
     }
     pMappingNode pNode = parentNode->next;
     parentNode->next = node;
     node->next = pNode;
     length++;
 }
 
 pMappingNode PQueue::next()
 {
     if (length==0)
         return NULL;
     pMappingNode oldHead = head;
     head = oldHead->next;
     length --;
     return oldHead;
 }
pMappingNode PQueue::traverse()
{
pMappingNode currentNode=head;
//cout<<"nodes in the queue are ";
while(currentNode!=NULL)
	{
//	cout<<"node-";
	int n=currentNode->stage;
//	for(int i=0;i<n;i++)
//	cout<< currentNode->mappingSequency[i]<<"->";
	currentNode=currentNode->next;
//	cout<<", ";
	}
}
void bnb_map()
 {
//cout<<"HEllo"<<endl;
//cout<<"calling initialize()"<<endl;
     initialize();
//cout<<"Now the processes are sort according to their communicaition requirement and rank is assignes to each process."<<endl;
sorting_process();
//cout<<"Now appMatrix is build on the basis of ranks"<<endl;
build_appMatrix();
//cout<<"HEllo"<<endl;
     float minCost = MAX_VALUE;
     float minulc = MAX_VALUE;
     pPQueue Q = new PQueue();

 system("rm -f output111/*");
     if (exist_locked_pe())
     {
         // this ruins the symmetric structure of the system completely.
         // although for some corner cases, symmetry still exists, we don't
         // consider it here.
         for (int i=0; i<numrows; i++)
         {
             for (int j=0; j<g_edge_size; j++)
             {
//cout<<"Creating the first level nodes"<<endl;
                 pMappingNode pNode = new MappingNode(i*g_edge_size+j);
                 if (pNode->is_illegal())
                     delete pNode;
                 else
                     Q->insert(pNode);
             }
         }
     }
     else
     {
         // To exploit the symmetric structure of the system, we only need
         // to map the first processes to one corner of the chip, as shown
         // in the following code.
         /*****************************************************************
          * And if we need to synthesize the routing table, then there is
          * not much symmetry property to be exploited
          *****************************************************************/
         
             int size = (g_edge_size+1)/2;
ofstream fout;
fout.open("outputnew11.txt", ios::trunc);
fout<<endl<<"Symmetric first level nodes are: ";
	    for(int snum=0;snum<slice_num;snum++)
	     {
	        for (int i=0; i<numrows; i++)
	       {
	           for (int j=0; j<g_edge_size; j++)
	           {
//fout<<size;
//cout<<"Creating the first level nodes";
int nodenum=i*g_edge_size+j+(s_size*snum) ;
fout<<"node-"<<nodenum<<" ";			
                     pMappingNode pNode = new MappingNode(i*g_edge_size+j+(s_size*snum));
                     if (pNode->is_illegal())
                         delete pNode;
                     else
                        Q->insert(pNode);

int len;
len=Q->Length();
//cout<<"Current Queue length is: "<< len<<endl;         
                 }
             }
          }
  fout.close();      
      }
 
     pMappingNode bestMapping = NULL;
     int min_ulc_hit_cnt = 0;
 
     while (!Q->empty())
     {
//counting++;
//if(counting<1)
//{
//cout<<endl<<"before calling  next() ";
Q->traverse();
//cout<<endl;
         pMappingNode pNode = Q->next();
//cout<<"after calling  next() ";
Q->traverse();
//cout<<endl;
int n=pNode->stage;
//cout<<"Node to be explored is node-";
//for(int i=0;i<n;i++)
//cout<< pNode->mappingSequency[i]<<"->";
//cout<<endl;
//cout<<"This node has llc = "<<pNode->llc<<" and ulc = "<<pNode->ulc<<endl<<pq_size;
         if (pNode->cost > minCost || pNode->llc > minulc)
         {
 //	cout<<"l";	
             delete pNode;
             continue;
         }

 //	cout<<"a";
         bool insertAllFlag = false;
//	 	cout<<"b";
         int prev_insert = 0;
//	cout<<"Current Queue length = "<<Q->Length()<<endl;
//	cout<<"Size of priority queue is: "<<pq_size;
	
/**********************************************************************
          *Change this to adjust the tradeoff between the solution quality     *
          *and the run time                                                    *
          **********************************************************************/

        if (Q->Length() < pq_size)
         {
insertAll:
//cout<<"Inside InsertAll:"<<endl;
//cout<<"InsertAllFlag's Status:"<<insertAllFlag<<endl;
             if (pNode->ulc == minulc && minulc < MAX_VALUE
                     && min_ulc_hit_cnt <= min_hit_threshold)
                 insertAllFlag = true;
//cout<<"InsertAllFlag's Status:"<<insertAllFlag<<endl;
             for (int i=prev_insert; i<gTileNum; i++)
             {
                 if (pNode->Expandable(i))
                 {
int l=pNode->stage;
//cout<<"Creating child Node-";
//for(int k=0;k<l;k++)
//cout<<pNode->mappingSequency[k]<<"->";
//cout<<i<<endl;
                     pMappingNode child = new MappingNode(*pNode, i);
//if(child->llc>minulc)
//cout<<"True1";
//else
//cout<<"false1";
//if(child->cost>minCost)
//cout<<"True2";
//else
//cout<<"false2";
//if(child->cost==minCost && bestMapping)
//cout<<"True3";
//else
//cout<<"false3";
//if(child->is_illegal())
//cout<<"True4";
//else
//cout<<"false4";
                     if (child->llc>minulc || child->cost>minCost
                             || (child->cost==minCost && bestMapping)
                             || child->is_illegal())
			
                         delete child;
			
                    else
                     {
                         if (child->ulc < minulc)
                         {
                             minulc = child->ulc;
                           //  cout <<endl<< "Current minimum cost upper bound is "
                                  //<< minulc << endl;
                             min_ulc_hit_cnt = 0;
 
                             //some new stuff here: we keep the mapping with min ulc
                             /*if (routing_table_synthesis)
                             {
                                 if (child->ulc<minCost)
                                 {
                                     if (bestMapping)
                                         delete bestMapping;
                                     bestMapping = new MappingNode(*child);
                                    minCost = child->ulc;
                                 }
                                 else if (child->ulc<minulc && !bestMapping)
                                     bestMapping = new MappingNode(*child);
                             }*/
                         }

                         if (child->stage==gProcNum)
                         {

                             minCost = child->cost;
                             if (child->stage<gProcNum)
                                 minCost = child->ulc;
                             if (minCost < minulc)
                                 minulc = minCost;
                             //cout<<"Current minimum cost is "<<minCost<<endl;
                             if (bestMapping)
                                 delete bestMapping;
                             bestMapping = child;
                         }
                         else
                         {
                             Q->insert(child);
int len;
len=Q->Length();
//cout<<"Current Queue length is: "<< len<<endl;         
                             if (Q->Length()>=pq_size && !insertAllFlag)
                             {
                                 prev_insert = i;
                                 goto selective_insert;
                             }
                         }
                     }
                 }
             }
            continue;
         }
         else
         {
 selective_insert:
//cout<<"Inside selective insert:"<<endl;
int diff=abs(pNode->ulc-minulc);
//cout<<diff<<endl;
            if ((abs(pNode->ulc-minulc)==0.01) && minulc<MAX_VALUE
                     && min_ulc_hit_cnt <= min_hit_threshold)
             {
                 min_ulc_hit_cnt++;
//cout<<"min_ulc_hit_cnt = "<<min_ulc_hit_cnt<<endl;
                 goto insertAll;
             }
             // In this case, we only select one child which has the
             // smallest partial cost. However, if the node is currently
             // the one with the minulc, then its child which
             // is generated by the corresponding minulc is
             // also generated
             int index = pNode->bestCostCandidate();
//cout<<" and the bestCostCandidate is :"<<index<<endl;
int l=pNode->stage;
//cout<<"Creating child Node-";
//for(int k=0;k<l;k++)
//cout<<pNode->mappingSequency[k]<<"->";
//cout<<index<<endl;
             pMappingNode child = new MappingNode(*pNode, index);
             if (child->llc>minulc || child->cost>minCost
                     || (child->cost==minCost&&bestMapping)
                     || child->is_illegal())
			{
				int s=child->stage;
//				cout<<"node-";
//				for(int k=0;k<s;k++)
//				cout<<child->mappingSequency[k]<<"->";
//				cout<<" is deleted"<<endl;
                 		delete child;
			}
             else
             {
                 if (child->ulc < minulc - 0.01)
                 {
int diffnew=minulc-0.01;
//cout<<diffnew<<endl;
                     // In this case, we should also insert other children
                     delete child;
                     insertAllFlag = true;
//cout<<" Go To InsertAll from SelectiveInsert"<<endl;
                     goto insertAll;
                 }
                 if (child->stage==gProcNum || child->llc==child->ulc)
                 {
                     minCost = child->cost;
//		     cout<<"Current MinCost is "<<minCost<<endl;
                     if (child->stage<gProcNum)
                         minCost = child->ulc;
                     if (minCost<minulc)
                         minulc = minCost;
  //                   cout<<"Current stage is "<<child->stage<<endl;
    //                 cout<<"Current minimum cost is "<<minCost<<endl;
                     if (bestMapping)
                         delete bestMapping;
                     bestMapping = child;
                 }
                 else
                     Q->insert(child);
int len;
len=Q->Length();
//cout<<"Current Queue length is: "<< len<<endl;         
             }
 
             if (pNode->ulc > minulc || pNode->ulc==MAX_VALUE)
             {   //cout<<"node-";
		 int r=pNode->stage;
		 //for(int k=0;k<r;k++)
		 //cout<<pNode->mappingSequency[k]<<"->";
		 //cout<<" is deleted"<<endl;
                 delete pNode;
                 continue;
             }
 int idx;
idx=pNode->bestulcCandidate();
//cout<<endl;
int m=pNode->stage;
//cout<<"The index value comes from calling bestCostCandidate function and the value of Index is ";
//for(int k=0;k<m;k++)
//cout<<pNode->mappingSequency[k]<<"->";
//cout<<index<<endl;
//cout<<"pNode->bestulcCandidate() = ";
//for(int k=0;k<m;k++)
//cout<<pNode->mappingSequency[k]<<"->";
//cout<<idx<<endl;
            if (index == pNode->bestulcCandidate())
             {
//		 cout<<endl<<"Before deleting the node, ";
		 Q->traverse();
//		 cout<<endl<<"node-";
		 int r=pNode->stage;
//		 for(int k=0;k<r;k++)
//		 cout<<pNode->mappingSequency[k]<<"->";
//		 cout<<" is deleted"<<endl;
                 delete pNode;
//		 cout<<"After deleting the node, ";
		 Q->traverse();
//cout<<"Control before continue";
                 continue;
             }
 
             index = pNode->bestulcCandidate();
#ifdef DEBUG
             if (!pNode->Expandable(index))
             {
                 cerr<<"Error in expanding at stage "<<pNode->Stage()<<endl;
                 cerr<<"index = "<<index<<endl;
                exit(-1);
             }
 #endif //DEBUG

//		cout<<"Creating child Node-";
                 int t;
		 t=pNode->stage;
//		for(int k=0;k<t;k++)
//		cout<<pNode->mappingSequency[k]<<"->";
//		cout<<index<<endl;
             child = new MappingNode(*pNode, index);
             if (child->llc > minulc || child->cost > minCost)
			{
			 int d;
			 d=pNode->stage;
//			for(int k=0;k<d;k++)
//			cout<<"node-"<<pNode->mappingSequency[k]<<"->";
//			cout<<index<<" is deleted"<<endl;
        	         delete child;
			}
             else
             {
                 if (child->ulc < minulc)
                 {
                     minulc = child->ulc;
  //                   cout<<"Current minimum cost upper bound is "<<minulc<<endl;
                     min_ulc_hit_cnt = 0;
                 }
                 if (child->stage==gProcNum || child->llc==child->ulc)
                 {
                     if (minCost==child->cost && bestMapping)
                         delete child;
                     else
                     {
                         minCost = child->cost;
	//		 cout<<"Current MinCost is "<<minCost<<endl;
                         if (child->stage<gProcNum)
                             minCost = child->ulc;
                         if (minCost<minulc)
                             minulc = minCost;
          //               cout<<"Current minimum cost is "<<minCost<<endl;
		//	 cout<<"Current stage is "<<child->stage<<endl;
                     
                         if (bestMapping)
                             delete bestMapping;
                         bestMapping = child;
                     }
                 }
                 else
			{
                     	Q->insert(child);
			int len=Q->Length();
		//	cout<<"Current Queue length is = "<<len<<endl;
			}

            }
         }
		//		cout<<"node-";
				int r=pNode->stage;
		//		for(int k=0;k<r;k++)
		//		cout<<pNode->mappingSequency[k];
		//		cout<<" is deleted"<<endl;
	
    
         delete pNode;
Q->traverse();
     //}
}
     cout<<"Totally "<<MappingNode::cnt<<" have been generated"<<endl;
    //int x= bestMapping->stage;
//for(int i=0;i<x;i++)
     if (bestMapping)
     {
         BBMMap(bestMapping);
         delete bestMapping;
     }
     else
         cout<<endl<<"Can not find a suitable solution."<<endl;
	
     BBMClear();
}
void BBMMap(pMappingNode bestMapping) {
rname=string(rem) + string(path) + string("application.config");
system(rname.c_str());
rname=string(rem) + string(path3) + string("*");
system(rname.c_str());
string applicationconfig =  string(path) + string("application.config"); 
string mapping_cost=string(path) + string("mapping_cost.txt"); 
//system("rm -f bnb/application.config");
//system("rm -f bnb/traffic/*");
ofstream config;
ofstream traffic;
ofstream mapping;
mapping.open(mapping_cost.c_str(), ios::trunc);
mapping.close();
char str_id[500];
char str_id1[500];
int size;
    for (int i=0; i<gProcNum; i++) 
        gProcess[i]->mapToTile(-1);
    //for (int i=0; i<gTileNum; i++) 
      //  gTile[i].mapToProc(-1);
    for (int i=0; i<gProcNum; i++) {
mapping.open(mapping_cost.c_str(), ios::app);
        int procId = app_rank_array[i];
        cout<<"Process "<<procId<<" is mapped to tile-"<<gProcess[procId]->mapToTile(bestMapping->mapToTile(i))<<endl;
mapping<<"Process "<<procId<<" is mapped to tile-"<<gProcess[procId]->mapToTile(bestMapping->mapToTile(i))<<endl;
mapping.close();
config<<"for Process "<<procId<<" : "<< t[procId].dst.size()<<endl;
if(t[procId].dst.size()==0)
continue;
config.open(applicationconfig.c_str(), ios::app);
config<<gProcess[procId]->mapToTile(bestMapping->mapToTile(i))<<" "<<"BWCBR"<<t[procId].dst.size()<<".so"<<endl;
config.close();
//        gTile[bestMapping->mapToTile(i)].mapToProc(procId);
    }
for(int i=0;i<gProcNum;i++)
{
int procId=app_rank_array[i];
int tile_id=gProcess[procId]->mapToTile(bestMapping->mapToTile(i));
sprintf(str_id, "%d", tile_id);
size=t[procId].dst.size();
for(int k=0;k<size;k++)
{
for(vector<int>::iterator m=t[procId].dst.begin();m!=t[procId].dst.end();m++)
{
sprintf(str_id1, "%d", k);
string traffic_filename = string(path3) + string(str_id) + string("_") + string(str_id1);
int z=*m;
int tid;
for(int s=0;s<gProcNum;s++)
{
if(app_rank_array[s]==z)
	{
	tid=gProcess[z]->mapToTile(bestMapping->mapToTile(s));
        break;
	}
else
continue;
}	
traffic.open(traffic_filename.c_str());
//traffic<<z<<endl;
traffic<<"PKT_SIZE"<<" "<<"8"<<endl;
traffic<<"DESTINATION FIXED"<<" "<<tid<<endl;
traffic<<"FLIT_ITNTERVAL"<<" "<<"2"<<endl;
traffic.close();
k++;
}
}
for(int k=0;k<size;k++)
{
for(vector<int>::iterator l=t[procId].load.begin();l!=t[procId].load.end();l++)
{
sprintf(str_id1, "%d", k);
string traffic_filename = string(path3) + string(str_id) + string("_") + string(str_id1);
traffic.open(traffic_filename.c_str(), ios::app);
traffic<<"LOAD"<<" "<<*l<<endl;
traffic.close();
k++;
}
}
}
int cost=0;
for(int i=0;i<gProcNum;i++)
{
for(int j=0;j<i;j++)
{
if(i==j)
continue;
int tile1=bestMapping->mappingSequency[i];
int tile2=bestMapping->mappingSequency[j];
cost+=appMatrix[i][j]*icnEnergyMatrix[tile1][tile2];
}
}
mapping.open(mapping_cost.c_str(), ios::app);
mapping<<"ulc is : "<<bestMapping->ulc<<endl;
mapping<<"Cost is : "<<cost<<endl;
cout<<endl;
cout<<"ulc is : "<<bestMapping->ulc<<endl;
cout<<"Cost is : "<<cost<<endl;
generate_rtable(bestMapping);
}
void BBMClear()
 {
     cout<<"Clear for Branch-and-Bound mapping"<<endl;
     if (icnEnergyMatrix)
     {
         for (int i=0; i<gTileNum; i++)
             delete []icnEnergyMatrix[i];
         delete []icnEnergyMatrix;
     }
     if (appMatrix)
     {
         for (int i=0; i<gProcNum; i++)
             delete []appMatrix[i];
         delete []appMatrix;
     }
     if (app_rank_array)
         delete []app_rank_array;
 }
void RANDOMMClear()
 {
     cout<<"Clear for Branch-and-Bound mapping"<<endl;
     if (icnEnergyMatrix)
     {
         for (int i=0; i<gTileNum; i++)
             delete []icnEnergyMatrix[i];
         delete []icnEnergyMatrix;
     }
     if (appMatrix)
     {
         for (int i=0; i<gProcNum; i++)
             delete []appMatrix[i];
         delete []appMatrix;
     }
     if (app_rank_array)
         delete []app_rank_array;
 }
bool exist_locked_pe()
{
return false;
} 
int MappingNode::bestCostCandidate()
 {
//cout<<"Inside bestCostCandidate()"<<endl;
     float minimal = MAX_VALUE;
     for (int i=0; i<gTileNum; i++)
         tileOccupancyTable[i] = false;
     for (int i=0; i<stage; i++)
{
         tileOccupancyTable[mappingSequency[i]] = true;
//cout<<"tileOccupancyTable["<<mappingSequency[i]<<"] = true"<<endl;
 }
     int index = -1;
     for (int tileId=0; tileId<gTileNum; tileId++)
     {
       if (tileOccupancyTable[tileId])
             continue;
         float additionalCost = 0;
         for (int i=0; i<stage; i++)
         {
             int tile1 = tileId;
             int tile2 = mappingSequency[i];
             additionalCost += appMatrix[i][stage] * icnEnergyMatrix[tile1][tile2];
//cout<<endl<<"with tile-"<<tileId<<" additional cost is"<<additionalCost<<endl;
             if (additionalCost >= minimal)
                 break;
         }
         if (additionalCost < minimal)
            minimal = additionalCost;
         index = tileId;
//cout<<"Current minimal cost is "<<minimal<<" and tiled is tile- "<<index<<endl;
     }
     return index;
 }
int MappingNode::bestulcCandidate()
 {
//cout<<"Inside bestulcCandidate() ";
     return mappingSequency[stage];
 }
bool MappingNode::Expandable(int tileId)
{
     //If it's an illegal mapping, then just return false
     for (int i=0; i<stage; i++)
     {
         if (mappingSequency[i] == tileId) //the tile has already been occupied
            return false;
     }
     return true;
 }
bool MappingNode::fixed_verify_BW_usage() {
    int *link_BW_usage_temp = new int[gLinkNum];
    memcpy(link_BW_usage_temp, link_BW_usage, sizeof(int)*gLinkNum);

    for (int i=0; i<stage; i++) {
        int tile1 = mappingSequency[i];
        int proc1 = app_rank_array[i];
        for (int j=stage; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            int proc2 = app_rank_array[j];
            if (gProcess[proc1]->outgoing_bw_requirement[proc2]) {
                for (unsigned int i=0; i<link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = link_usage_list[tile1][tile2][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->outgoing_bw_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > link_bw[link_id]) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }

            if (gProcess[proc1]->incoming_bw_requirement[proc2]) {
                for (unsigned int i=0; i<link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = link_usage_list[tile2][tile1][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->incoming_bw_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > link_bw[link_id]) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }
        }
    }
    for (int i=stage; i<gProcNum; i++) {
        int tile1 = mappingSequency[i];
        int proc1 = app_rank_array[i];
        for (int j=i+1; j<gProcNum; j++) {
            int tile2 = mappingSequency[j];
            int proc2 = app_rank_array[j];
            if (gProcess[proc1]->outgoing_bw_requirement[proc2]) {
                for (unsigned int i=0; i<link_usage_list[tile1][tile2].size(); i++) {
                    int link_id = link_usage_list[tile1][tile2][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->outgoing_bw_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > link_bw[link_id]) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }

            if (gProcess[proc1]->incoming_bw_requirement[proc2]) {
                for (unsigned int i=0; i<link_usage_list[tile2][tile1].size(); i++) {
                    int link_id = link_usage_list[tile2][tile1][i];
                    link_BW_usage_temp[link_id] += gProcess[proc1]->incoming_bw_requirement[proc2];
                    if (link_BW_usage_temp[link_id] > link_bw[link_id]) {
                        delete []link_BW_usage_temp;
                        return false;
                    }
                }
            }
        }
    }
    delete []link_BW_usage_temp;
    return true;
}
bool parse_apcg(char * fileName) {
int max=0;
int load;
 
    char inputLine[MAX_LINE];
    ifstream inputFile(fileName);
    int src, dst, commVol, BW;
    if (!inputFile) {
        cerr<<"Error in openning file "<<fileName<<endl;
        return false;
    }
	 
    while (inputFile.getline(inputLine, MAX_LINE, '\n')) {
	
        if (strchr(inputLine, '#'))
            continue;
        if (strlen(inputLine) == 0)
            continue;
        if (inputLine[0] == '\r' || inputLine[0] == '\n')
            continue;
        if (!sscanf(inputLine, "%d\t%d\t%d\t%d", &src, &dst, &commVol, &BW)) 
            return false;
    if(commVol>max)
	max=commVol;
	t[src].dst.push_back(dst);
        t[src].toVolume.push_back(commVol);
         gProcess[src]->toComm[dst] = commVol;
         gProcess[src]->outgoing_bw_requirement[dst] = BW;
         gProcess[dst]->fromComm[src] = commVol;
         gProcess[dst]->incoming_bw_requirement[src] = BW; 
}
for(int i=0;i<gProcNum;i++)
{
if(t[i].dst.size()==0)
continue;
for(vector<int>::iterator j=t[i].toVolume.begin(); j!=t[i].toVolume.end();j++)
{
if(load_condition=='f')
load = 100;
else if(load_condition=='o')
load=((*j)*100)/max;
t[i].load.push_back(load);
}
}
inputFile.close();
 return true;
}
void sorting_process() {
    //sort them according to the sum of each process's ingress and egress 
    //communication volume
    for (unsigned int i=0; i<gProcess.size(); i++) {
        gProcess[i]->totalCommVol = 0;
        for (unsigned int k=0; k<gProcess.size(); k++) {
            gProcess[i]->totalCommVol += gProcess[i]->toComm[k];
            gProcess[i]->totalCommVol += gProcess[i]->fromComm[k];
        }
    }
    // Now rank them
    int cur_rank = 0;
    // locked PEs have the highest priority
    for (unsigned int i=0; i<gProcess.size(); i++) {
        /*if (gProcess[i]->is_locked()) {
            app_rank_array[cur_rank] = i;
            gProcess[i]->rank = cur_rank ++;
        }
        else*/
            gProcess[i]->rank = -1;
    }
    // the remaining PEs are sorted based on their comm volume
    for (unsigned int i=cur_rank; i<gProcess.size(); i++) {
        int max = -1;
        int maxid = -1;
        for (int k=0; k<gProcNum; k++) {
            if (gProcess[k]->rank != -1)
                continue;
            if (gProcess[k]->totalCommVol > max) {
                max = gProcess[k]->totalCommVol;
                maxid = k;
            }
        }
        gProcess[maxid]->rank = i;
        app_rank_array[i] = maxid;
    }
for(int i=0;i<gProcNum;i++)
cout<<"app_rank_array["<<i<<"] = "<<app_rank_array[i]<<endl;
}
void build_appMatrix() {
    appMatrix = new int*[gProcNum];
    for (int i=0; i<gProcNum; i++) 
        appMatrix[i] = new int[gProcNum];
    //fill it with corresponding value
    for (int i=0; i<gProcNum; i++) {
        int row = gProcess[i]->rank;
        for (int j=0; j<gProcNum; j++) {
            int col = gProcess[j]->rank;
            appMatrix[row][col] = gProcess[i]->fromComm[j] + gProcess[i]->toComm[j];
        }
    }
    //Sanity checking
#ifdef DEBUG
    for (int i=0; i<gProcNum; i++) {
        for (int j=0; j<gProcNum; j++) {
            if (appMatrix[i][j]<0) {
                cerr<<"Error for <0"<<endl;
                exit(1);
            }
            if (appMatrix[i][j]!=appMatrix[j][i]) {
                cerr<<"Error. The process matrix is not symetric."<<endl;
                exit(1);
            }
        }
    }
#endif //DEBUG
for (int i=0; i<gProcNum; i++)
        for (int j=0; j<gProcNum; j++) 
            cout<<"appMatrix["<<i<<"]["<<j<<"] = "<<appMatrix[i][j]<<endl;
}
void build_icnEnergyMatrix() {
int linkId;
    icnEnergyMatrix = new float*[gTileNum];
    for (int i=0; i<gTileNum; i++) 
        icnEnergyMatrix[i] = new float[gTileNum];
    for (int src=0; src<gTileNum; src++) {
        for (int dst=0; dst<gTileNum; dst++) {

cout<<"Path from tile-"<<src<<" to tile-"<<dst<<" consists of links- ";
            float energy = 0;
            pTile currentTile = &gTile[src];
            energy += currentTile->Cost();
if(src==dst)
{
linkId=-1;
cout<<linkId;
}
else{
            while (currentTile->GetId()!= dst) {
                linkId = currentTile->RouteToLink(src, dst);
                cout<<linkId<<" ";
                pLink pL = &gLink[linkId];
                energy += pL->Cost();
                currentTile = &gTile[pL->ToTile()];
                //cout<<currentTile->GetId();
                energy += currentTile->Cost();
            }
}
cout<<endl;
            icnEnergyMatrix[src][dst] = energy;
        }
    }
}
void generate_rtable(pMappingNode bestMapping)
{
string rname=string(rem) + string(path1) + string("*");
system(rname.c_str());
create_files();
ofstream rtable;
ifstream rtable1;
int size,toid;
int pathid=0;
char str_id[500];
char str_id1[500];
for(int i=0;i<gProcNum;i++)
{
int count=0;
int procId=app_rank_array[i];
int tile_id=gProcess[procId]->mapToTile(bestMapping->mapToTile(i));
sprintf(str_id, "%d", tile_id);
size=t[procId].dst.size();
if(size==0)
continue;
string rtable_filename = string(path1) + string(str_id) + string(".txt");
for(vector<int>::iterator m=t[procId].dst.begin();m!=t[procId].dst.end();m++)
{
int y=*m;
for(int s=0;s<gProcNum;s++)
{
if(app_rank_array[s]==y)
	{
	toid=gProcess[y]->mapToTile(bestMapping->mapToTile(s));
        break;
	}
else
continue;
}
int linkId,tilid;
		pTile nextTile = &gTile[tile_id];
                while (nextTile->GetId()!= toid) 
		{
tilid=nextTile->GetId();
                linkId = nextTile->RouteToLink(tile_id, toid);         
                pLink pL = &gLink[linkId];
		nextTile = &gTile[pL->ToTile()];
sprintf(str_id1, "%d", tilid);
string rtable_filename1 = string(path1) + string(str_id1) + string(".txt");
rtable.open(rtable_filename1.c_str(), ios::app);
rtable<<tile_id<<" "<<toid<<" "<<nextTile->GetId()<<" "<<pathid<<endl;
rtable.close();
}
}
}
remove_unwanted_files();
}
void generate_topology_IR_config()
{
ofstream topology;
rname=string(rem) + string(path2) + string("topology_TD.config");
system(rname.c_str());
string topologyconfig=string(path2) + string("topology_TD.config");
topology.open(topologyconfig.c_str(), ios::app);
int a;
Link *glink;
for(int i=0;i<gTileNum;i++)
{
topology<<gTile[i].id<<" ";
for(int j=0;j<gTile[i].GetGoLinkNum();j++)
	{
          a=gTile[i].GoLink(j)->ToTile();
//a=glink->ToTile();
//a=(*glink).ToTile();
topology<<a<<" ";
	}
topology<<"-1";
topology<<endl;
}
topology.close();
}

void generate_link_length()
{
rname=string(rem) + string(path2) + string("link_length");
system(rname.c_str());
string linklength=string(path2) + string("link_length");
ofstream file;
int a;
int * b= NULL;
int size;
Link *glink;
int length;
int count;
for(int i=0;i<gTileNum;i++)
{
   size=gTile[i].GetGoLinkNum();
   b=new int[size];
	for(int j=0;j<gTile[i].GetGoLinkNum();j++)
	{
          a=gTile[i].GoLink(j)->ToTile();
	  b[j]=a;
	}
	for(int k=0;k<gTileNum;k++)
	{
	count=0;
	for(int l=0;l<gTile[i].GetGoLinkNum();l++)
	{
	if(k==b[l])
	{
	length=link_length;
	file.open(linklength.c_str(), ios::app);
	file<<i<<" "<<k<<" "<<length<<endl;
	file.close();
	count++;
	break;
	}
	}
	if(count==0)
	{
	length=-1;
	file.open(linklength.c_str(), ios::app);
	file<<i<<" "<<k<<" "<<length<<endl;
	file.close();
	}
	}
}
}
void create_files()
{
ofstream file;
char str_id1[500];
for(int i=0;i<gTileNum;i++)
{
sprintf(str_id1, "%d", i);
string rtable_filename1 = string(path1) + string(str_id1) + string(".txt");
file.open(rtable_filename1.c_str(), ios::app);
file<<"NUM_RTABLE_ENTRIES"<<" "<<" "<<" "<<" "<<" "<<" "<<" "<<endl;
file<<"#"<<"Source"<<" "<<"destination"<<" "<<"nexttile"<<" "<<"pathid"<<endl;
file.close();
}
}
void remove_unwanted_files()
{
int number_of_lines=0, entry;
ifstream fp;
string s;
ofstream file;
char str_id1[500];
for(int i=0;i<gTileNum;i++)
{
number_of_lines=0;
sprintf(str_id1, "%d", i);
string rtable_filename1 = string(path1) + string(str_id1) + string(".txt");
fp.open(rtable_filename1.c_str(), ios::in);
while (getline(fp, s))
{
number_of_lines++;
}
fp.close();
cout<<number_of_lines<<endl;
entry=number_of_lines-2;
if(number_of_lines==2)
{
//remove(rtable_filename1.c_str());
file.open(rtable_filename1.c_str(), ios::in | ios::out);
file.seekp(19, ios::beg);
file<<"0";
file.close();
}
else
{
file.open(rtable_filename1.c_str(), ios::in | ios::out);
file.seekp(19, ios::beg);
file<<entry;
file.close();
}
}
}
/////////////////////////////Random_Mapping/////////////////////////////
void randommapping()
{
pMappingNode randomMapping=NULL;
initialize();
int i=0;
for(int i=0;i<gProcNum;i++)
{
app_rank_array[i]=i;
gProcess[i]->rank=i;
}
build_appMatrix();
pMappingNode pNode = new MappingNode(i);
for(int i=0;i<gTileNum;i++)
pNode->mappingSequency[i]=i;
randomMapping=pNode;
randommap(randomMapping);
}
/////////////////////////////Random_Mapping_1//////////////////////////////////////
void randommapping1()
{
pMappingNode randomMapping=NULL;
initialize();
int nodenum=0;
int n,a,count=0;
srand(time(0));
int num[1000];
for(int i=0;i<gProcNum;i++)
{
app_rank_array[i]=i;
gProcess[i]->rank=i;
}
build_appMatrix();
pMappingNode pNode = new MappingNode(nodenum);
for(int i=0;i<gTileNum;)
{
a=i;
if(i==0)
count++;
n=rand() % gTileNum ;
for(int j=0;j<i;j++)
{
if(n==num[j])
{
count=0;
break;
}
else
{
count++;
}
}
if(count>0)
{
num[i]=n;
i++;
count=0;
}
else 
i=a;
}
for(int k=0;k<gTileNum;k++)
pNode->mappingSequency[k]=num[k];
randomMapping=pNode;
randommap(randomMapping);
}
void randommap(pMappingNode randomMapping) {
rname=string(rem) + string(path) + string("application.config");
system(rname.c_str());
rname=string(rem) + string(path3) + string("*");
system(rname.c_str());
string applicationconfig =  string(path) + string("application.config"); 
string mapping_cost=string(path) + string("mapping_cost.txt"); 
ofstream config;
ofstream traffic;
ofstream mapping;
mapping.open(mapping_cost.c_str(), ios::trunc);
mapping.close();
char str_id[500];
char str_id1[500];
int size;
    for (int i=0; i<gProcNum; i++) 
        gProcess[i]->mapToTile(-1);
        for (int i=0; i<gProcNum; i++) 
 {
        int procId = app_rank_array[i];
mapping.open(mapping_cost.c_str(), ios::app);
        cout<<"Process "<<procId<<" is mapped to tile-"<<gProcess[procId]->mapToTile(randomMapping->mapToTile(i))<<endl;
mapping<<"Process "<<procId<<" is mapped to tile-"<<gProcess[procId]->mapToTile(randomMapping->mapToTile(i))<<endl;
mapping.close();
config<<"for Process "<<procId<<" : "<< t[procId].dst.size()<<endl;
if(t[procId].dst.size()==0)
continue;
config.open(applicationconfig.c_str(), ios::app);
config<<gProcess[procId]->mapToTile(randomMapping->mapToTile(i))<<" "<<"BWCBR"<<t[procId].dst.size()<<".so"<<endl;
config.close();
//        gTile[randomMapping->mapToTile(i)].mapToProc(procId);
    }
for(int i=0;i<gProcNum;i++)
{
int procId=app_rank_array[i];
int tile_id=gProcess[procId]->mapToTile(randomMapping->mapToTile(i));
sprintf(str_id, "%d", tile_id);
size=t[procId].dst.size();
for(int k=0;k<size;k++)
{
for(vector<int>::iterator m=t[procId].dst.begin();m!=t[procId].dst.end();m++)
{
sprintf(str_id1, "%d", k);
string traffic_filename = string(path3) + string(str_id) + string("_") + string(str_id1);
int z=*m;
int tid;
for(int s=0;s<gProcNum;s++)
{
if(app_rank_array[s]==z)
	{
	tid=gProcess[z]->mapToTile(randomMapping->mapToTile(s));
        break;
	}
else
continue;
}	
traffic.open(traffic_filename.c_str());
//traffic<<z<<endl;
traffic<<"PKT_SIZE"<<" "<<"8"<<endl;
traffic<<"DESTINATION FIXED"<<" "<<tid<<endl;
traffic<<"FLIT_ITNTERVAL"<<" "<<"2"<<endl;
traffic.close();
k++;
}
}
for(int k=0;k<size;k++)
{
for(vector<int>::iterator l=t[procId].load.begin();l!=t[procId].load.end();l++)
{
sprintf(str_id1, "%d", k);
string traffic_filename = string(path3) + string(str_id) + string("_") + string(str_id1);
traffic.open(traffic_filename.c_str(), ios::app);
traffic<<"LOAD"<<" "<<*l<<endl;
traffic.close();
k++;
}
}
}
int cost=0;
for(int i=0;i<gProcNum;i++)
{
for(int j=0;j<i;j++)
{
if(i==j)
continue;
int tile1=randomMapping->mappingSequency[i];
int tile2=randomMapping->mappingSequency[j];
cost+=appMatrix[i][j]*icnEnergyMatrix[tile1][tile2];
}
}
mapping.open(mapping_cost.c_str(), ios::app);
mapping<<"Cost is : "<<cost<<endl;
mapping.close();
cout<<endl;
//cout<<"ulc is : "<<randomMapping->ulc<<endl;
cout<<"Cost is : "<<cost<<endl;
generate_rtable(randomMapping);
RANDOMMClear();
}
main()
{
cout<<"hello inside main"<<endl;
cout<<"Enter number of rows : ";
cin>>numrows;
cout<<endl;
cout<<"Enter number of columns : ";
cin>>g_edge_size;
cout<<endl;
cout<<"Enter number of slices : ";
cin>>slice_num;
cout<<endl;
cout<<"Enter number of processes to be mapped : ";
cin>>gProcNum;
cout<<endl;
cout<<"Enter o for Optimized mapping and r for Random Mapping : ";
cin>>type;
cout<<endl;
snum=slice_num;
if(type=='r')
{
cout<<"Note : Two vesrions are available for random mapping...."<<endl;
cout<<"Type 1 for version 1 and 2 for Version 2"<<endl;
cin>>ver;
cout<<endl;
}
cout<<"Enter name of the file having application characterisitcs : ";
cin>>apcgfilename;
cout<<endl;
cout<<"Specify type of load you want by pressing o for optimized load and f for full load : ";
cin>>load_condition;
cout<<endl;
dimension[0]=numrows;
dimension[1]=g_edge_size;
dimension[2]=slice_num;
if(type=='o')
folder_name=string("bnb_results_extra_111/optimized-");
else if(type=='r')
folder_name=string("bnb_results_extra_111/random-");
rtable=string("rtable");
topology=string("topology");
traff=string("traffic");
for(int i=0;i<=2;i++)
{
sprintf(folder,"%d",dimension[i]);
if(i<2)
folder_name+=string(folder) + string("*");
else
folder_name+=string(folder);
}
mkdir(folder_name.c_str(), S_IRWXU);
path=string(folder_name)+string("/");
path1=string(path) + string(rtable) + string("/");
path2=string(path) + string(topology) + string("/");
path3=string(path) + string(traff) + string("/");
mkdir(path1.c_str(), S_IRWXU);
mkdir(path2.c_str(), S_IRWXU);
mkdir(path3.c_str(), S_IRWXU);
s_size=numrows * g_edge_size; 
gTileNum=s_size*slice_num;
cout<<gTileNum<<endl;
if(gTileNum<gProcNum)
{
cout<<"Mapping can't be done!!!!!!!!!"<<endl;
cout<<"Note : Number of process can't be set more than the total number of tiles in architecture"<<endl;
exit(1);
}
int single=((2 * (g_edge_size-1) * numrows) +((numrows-1) * 2 * g_edge_size))*slice_num;
int between=(2*s_size)*(slice_num-1);
gLinkNum=single+between;
gTile = new Tile[gTileNum]();
gLink = new Link[gLinkNum]();
for (int i=0; i<gTileNum; i++)
{
cout<<"Tile-"<<i<<" has outgoing and incoming links : ";
gTile[i].AttachLink(gLink);
cout<<endl;
gTile[i].initialize_router();
}
if(type=='o')
bnb_map();
else if(type=='r')
{
if(ver==1)
randommapping();
else if(ver==2)
randommapping1();
}
generate_topology_IR_config();
generate_link_length();
return 0;
}
