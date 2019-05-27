/************************************************************************/
/* 

//.ctrGraph format

<file> := <#vertices>(<vertex>*)<#planes>(<plane>*)
<plane>:= 'p' <#plane param> #edges_number_on_plane_int'\n'(<edge>*)
<vertex> := 'v' vx_float vy_float vz_float
<edge> := 'e' v_ind1_int v_ind2_int m_ind1_int m_ind2_int
<comment> := '#' <character>*

before <#vertices> and <#planes> 'n' should be added.this is easier 
from programming view, so that whenever one line is scanner, first 
character is read and check if it's comment.

This format needs a well preprocessed contour. If this kind of file
is taken in, it will assume that consistency has already been handled
outside the project. That means, all intersection contour vertices have 
already been computer, and all the material configurations on different
planes are consistent too.

INDEX starts with 0, for nodes and edges too.
*/
/************************************************************************/

#include "../ContourHandler/ContourHandler.h"
#include "../math/mymath.h"
#include <string>
#include <set>
using namespace std;

vector<float> graph_nodes;
vector<float> graph_plane_params;
vector<vector<int> > graph_plane_edges;

/*--	read contours graphs in		--*/
void errorFileFormat()
{
	cout<<"Invalid file!"<<endl;
	cout<<"File format should be:"<<endl;
	cout<<"n #vergices\n"
		<<"v x1 y1 z1\n"
		<<"v x2 y2 z2\n"
		<<"...\n"
		<<"v xn yn zn\n"
		<<"n #planes\n"
		<<"p pa_1 pb_1 pc_1 pd_1 #edges_1"
		<<"e v1_1_ind v1_2_ind m1_1_ind m1_2_ind\n"
		<<"e v2_1_ind v2_2_ind m2_1_ind m2_2_ind\n"
		<<"...\n"
		<<"e vn1_1_ind vn1_2_ind mn1_1_ind mn1_2_ind\n"
		<<"p pa_2 pb_2 pc_2 pd_2 #edges_2"
		<<"...\n";
	cout<<"comments are allowed as a separate line starting with # , they are not allowed inside a line	or at the end of a line!\n";
}

void clearGraph()
{
	//clear the whole graph
	graph_nodes.clear();
	graph_plane_params.clear();
	for (unsigned int i = 0; i < graph_plane_edges.size(); i++)
		graph_plane_edges[ i ].clear();
	graph_plane_edges.clear();
}

bool readUntilLine_Prefix( ifstream& fin, char prefix )
{
    char chr;
	string str;

	if( fin.eof())
	{
		errorFileFormat();
		clearGraph();
		fin.close();
		return false;
	}

	do
	{	
		fin >> chr;
		if( chr == prefix )	//found!!!
			return true;
        getline( fin, str );	//skip the whole line
	}
	while( !fin.eof() );

	//coming here, outside of the fiel, still not able to return
	errorFileFormat();
	clearGraph();
	fin.close();
	return false;
}



void ContourHandler::readOneContour3_v3(char* filename, float bbox[ 6 ], bool& bboxset)
{
	ifstream fin ( filename );
	if( fin.bad())
	{
		cout<<"Unable to open file "<<filename<<" to read!"<<endl;
		return;
	}

	//read in nodes, and planes		
	int nodeNum;
	int planeNum;
	int edgeNum;

	//read until the #vertices is read in
	bool suc;
	suc = readUntilLine_Prefix(fin, 'n');
	if( !suc )
		return;
	fin >> nodeNum;

	//////////////////////////////////////////////////////////////////////////
	//cout<<"nodeNum:"<<nodeNum<<endl;
	//////////////////////////////////////////////////////////////////////////

	//read all the vertices in
	graph_nodes.resize( 3 * nodeNum);
	int pos = 0;
	for( int i = 0; i < nodeNum; i ++ )
	{
		suc = readUntilLine_Prefix( fin, 'v');
		if( !suc )	return;
		for( int j= 0; j < 3; j ++ )
			fin >> graph_nodes[ pos ++ ] ;
	}

	//////////////////////////////////////////////////////////////////////////
	/*cout<<"positions:\n";
	for( int i = 0; i <  nodeNum; i ++ )
		cout<<graph_nodes[ i * 3  ]<<","<<graph_nodes[ 3 * i +  1]<<","
		<<graph_nodes[ 3 * i + 2 ]<<" ";
	cout<<"\n===\n";*/
	//////////////////////////////////////////////////////////////////////////

	//read until 'n' is found, then read in #planes
	suc = readUntilLine_Prefix( fin, 'n' );
	if( !suc ) return;
	fin >> planeNum;

	//////////////////////////////////////////////////////////////////////////
	//cout<<"planenum:"<<planeNum<<endl;
	//////////////////////////////////////////////////////////////////////////
	
	//read all the planes in
	graph_plane_params.resize( 4 * planeNum );
	graph_plane_edges.resize( planeNum );
	pos = 0;
	int pos1 = 0;
	for( int i = 0; i < planeNum; i ++ )
	{
		//read until the 'p' line is found
		suc = readUntilLine_Prefix( fin, 'p' );
		if( !suc ) return;

		//read in plane parameters and edge number on it
		for( int j = 0; j < 4; j ++ )
			fin >> graph_plane_params[ pos ++ ];			
		fin >> edgeNum;
		
		//////////////////////////////////////////////////////////////////////////
		/*cout<<"4 parameters:"
			<<graph_plane_params[ pos - 4 ]<<" "
			<<graph_plane_params[ pos - 3 ]<<" "
			<<graph_plane_params[ pos - 2 ]<<" "
			<<graph_plane_params[ pos - 1 ]<<" "<<endl;*/
		//////////////////////////////////////////////////////////////////////////

		graph_plane_edges[ i ].resize( edgeNum * 4 );
		pos1 = 0;
		for( int j = 0; j < edgeNum; j ++ )
		{
			//read edge and material configuration in
            suc = readUntilLine_Prefix( fin, 'e' );
			if( !suc ) return;		

			for( int j=0; j < 4; j ++ )
				fin >> graph_plane_edges[ i ][ pos1 ++ ];				
			//////////////////////////////////////////////////////////////////////////
			/*	cout<<"edge:"
					<<graph_plane_edges[ i ][ pos1 - 4  ]<<" "
					<<graph_plane_edges[ i ][ pos1 - 3 ]<<" "
					<<graph_plane_edges[ i ][ pos1 - 2 ]<<" "
					<<graph_plane_edges[ i ][ pos1 - 1 ]<<endl;*/
			//////////////////////////////////////////////////////////////////////////

		}		
	}
	
	//set bounding box
	bbox[ 0 ] = bbox[ 3 ] = graph_nodes[ 0 ];
	bbox[ 1 ] = bbox[ 4 ] = graph_nodes[ 1 ];
	bbox[ 2 ] = bbox[ 5 ] = graph_nodes[ 2 ];
	bboxset = true;
	
	pos = 3;
	for( int i = 1; i < nodeNum; i ++ )
	{

		for( int j = 0; j < 3; j++ )
		{
			if( bbox[ j ] > graph_nodes[ pos ])	//bbox 0, 1, 2 saves the min in direction j. 3, 4,5 saves max
			{
				bbox[ j ] = graph_nodes[ pos ];
			}
			else if( bbox[ 3 + j ] < graph_nodes[ pos ])
			{
				bbox[ 3 + j ] = graph_nodes[ pos ];
			}
			pos ++;
		}		
	}	
}


/*		--	process the contour		--
1. ver2planelistpos, and ver2planelist
2.
 a. resave them in the conventional format, each plane has its own vertices, and contour edges
 b. for those planes share common line, save them as sets, for partition usage
*/

vector<vector<int> > ExistGraphNodeInPlane;	//for a node i, if it exists in plane j
vector<vector<int> > NodeIndex_Graph2Plane;	//for a node index i in graph, its index in plane j
vector<vector<int> > NodeIndex_Plane2Graph;	//for a node index i in plane i, its index in graph

void markExist()
{
	//////////////////////////////////////////////////////////////////////////
	/*cout<<"all the edges:\n";
	for( int i = 0; i < graph_plane_edges.size(); i ++ )
	{
		for( int j = 0; j < graph_plane_edges[ i ].size(); j++ )
		{
			cout<<graph_plane_edges[ i ][ j ]<<" ";
		}
		cout<<"\n-----\n";
	}
	cout<<endl;*/

	//////////////////////////////////////////////////////////////////////////

	int nodeNum = graph_nodes.size() / 3;
	
	ExistGraphNodeInPlane.resize( nodeNum );
	NodeIndex_Graph2Plane.resize( nodeNum );
	int planeNum = graph_plane_edges.size();
	for( int i = 0; i < nodeNum; i++ )
	{
		ExistGraphNodeInPlane [i ].resize( planeNum , 0);		
		NodeIndex_Graph2Plane[ i ].resize( planeNum, -1 );
	}
	
	//step 1. mark existence based on edge composition
	for( int i = 0; i < planeNum; i ++ )
	{
		int edgeNum = graph_plane_edges[ i ].size() / 4;
		int pos = 0;
		for( int j = 0; j < edgeNum ; j ++ )
		{
            ExistGraphNodeInPlane[ graph_plane_edges[ i ][ pos ] ][ i ] = 1;
			ExistGraphNodeInPlane[ graph_plane_edges[ i ][ pos + 1] ][ i ] = 1;
			pos += 4;
		}
	}
    
	//step 2. resave all the vertcies in each plane, set the map at the same time	
	NodeIndex_Plane2Graph.resize( planeNum );
	//planeCtrVer.resize( planeNum );
	for( int i = 0; i < planeNum; i ++ )
	{
		NodeIndex_Plane2Graph[ i ].reserve( nodeNum / planeNum );
		//push those vertices in
		int indInPlane = 0;
		for( int j = 0; j < nodeNum; j ++ )
		{
			if( ExistGraphNodeInPlane[ j ][ i ])	//node j exists in plane i
			{
				//plane 2 graph
				NodeIndex_Plane2Graph[ i ].push_back( j );	//j'th node
				//graph 2 plane
				NodeIndex_Graph2Plane[ j ][ i ] = indInPlane;	//the j'th node in plane i is at postion indInPlane
				indInPlane ++;
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////
	//cout<<"plane2graph node index:"<<endl;
	//for( int i = 0; i < planeNum; i ++ )
	//{
	//	//output all the indices on the planes
	//	for( int j = 0; j < NodeIndex_Plane2Graph[ i ].size(); j++  )
	//	{
	//		cout<<NodeIndex_Plane2Graph[ i ][ j ]<<" ";
	//	}
	//	cout<<"\n------\n";
	//}
	//////////////////////////////////////////////////////////////////////////

}

/*
according to the mapping information, find all the common vertices between two planes
sort all the common vertices out, merge those planes with the same set of 
intersection points, find the set of planes that intersect at the same set of points

*/
vector<vector<int> > interPtIndVec;
vector<vector<int> > radialPlanePairs;


//return -1, if first one is smaller, 0 if equal, and 1 if larger
int compare2Vector( vector<int>& vec1, vector<int>& vec2 )
{
	int num = vec1.size() < vec2.size() ? vec1.size() : vec2.size();	//minnum
	int maxnum = vec1.size() + vec2.size() - num;	//maxnum

	int pos = 0;
	while( pos < num && vec1[ pos ] == vec2[ pos ] )	//find the first that's different
		pos ++;

	int rslt;	
	if( pos != num )	//find one that's different
		rslt = ( vec1[ pos ] < vec2[ pos ] ) ? -1 : 1;
	else	//pos == num, equal until num
		rslt = ( num == maxnum ) ? 0 : ( (vec1.size() < vec2.size())? -1 : 1);
	return rslt;
}

void radialSetAnalyze()
{
	//step1. find all the intersection points between every two planes
	int planeNum = graph_plane_edges.size();
    int size = planeNum * ( planeNum - 1)/2;
	
	interPtIndVec.resize( size );
	radialPlanePairs.resize( size );

	int pos = -1;
	for( int i = 0; i < planeNum - 1; i ++ )
	{
		for( int j = i + 1; j < planeNum ; j ++ )
		{
			pos ++;
			radialPlanePairs[ pos ].resize(2 );
			radialPlanePairs[ pos ][ 0 ] =  i ;
			radialPlanePairs[ pos ][ 1 ] =  j ;

			interPtIndVec.reserve( 2 * planeNum );
			//go through all the graph nodes on plane i and j, if same, push into interPtIndSet
			int ptr1, ptr2;
			ptr1 = ptr2 = 0;
			int node1, node2;
			node1 = NodeIndex_Plane2Graph[ i ].size();
			node2 = NodeIndex_Plane2Graph[ j ].size();

			//this algorithm takes advantage that all the crpding node index is
			//in increasing order in the plane2graph node map.
			while( ptr1 < node1 && ptr2 < node2 )
			{
				//move the first pointer, until the same index is found
				if( NodeIndex_Plane2Graph[ i ][ ptr1 ] < NodeIndex_Plane2Graph[ j ][ ptr2 ])
				{
					while( ( NodeIndex_Plane2Graph[ i ][ ptr1 ] < NodeIndex_Plane2Graph[ j ][ ptr2 ] )
						&& ( ptr1 < node1 )	)	//move the first pointer until the value is >= the one in second one.
						ptr1 ++;								
				}
				//move the second pointer, until the same index is found
				else if( NodeIndex_Plane2Graph[ i ][ ptr1 ] > NodeIndex_Plane2Graph[ j ][ ptr2 ] )
				{
					while( ( NodeIndex_Plane2Graph[ i ][ ptr1 ] > NodeIndex_Plane2Graph[ j ][ ptr2 ] )
						&& ( ptr2 < node2 )	)	//move the first pointer until the value is >= the one in second one.
						ptr2 ++;
				}

				if( ptr1 >= node1 || ptr2 >= node2 )
					break;

				//if after moving, push back the same index
				if(NodeIndex_Plane2Graph[ i ][ ptr1 ] == NodeIndex_Plane2Graph[ j ][ ptr2 ] )
				{
					//the two are the same, push back the nodes into the interPt
					interPtIndVec[ pos ].push_back( NodeIndex_Plane2Graph[ i ][ ptr1 ] );
					ptr1 ++;
					ptr2 ++;
				}		
			}            
		}
	}

	//step 2, sort the two vectors based on the intersection points indices
	//take advantage of the fact that the intersection points index set are all in increasing orders.
	//sort them in alphabetical order	: select sort
    for( int i = 0; i < size - 1; i++ )
	{
		int minInd = i;

		for( int j = i + 1; j < size; j ++ )
		{
			//check if current vector is smaller than the one at minInd
			if( compare2Vector( interPtIndVec[ minInd ], interPtIndVec[ j ]) == 1) //the first one is larger
				minInd = j;
		}

		if( minInd == i )
			continue;

		//swap the one at minInd and i
		vector<int> temp( interPtIndVec[ i ].size() );
		for (unsigned int j = 0; j < temp.size(); j++)
			temp[ j ] = interPtIndVec[ i ][ j ];
		interPtIndVec[ i ].resize( interPtIndVec[ minInd].size() );
		for (unsigned int j = 0; j < interPtIndVec[minInd].size(); j++)
			interPtIndVec[ i ][ j ] = interPtIndVec[ minInd ][ j ];
		interPtIndVec[ minInd  ].resize( temp.size() );
		for (unsigned int j = 0; j < temp.size(); j++)
			interPtIndVec[ minInd ][ j ] = temp[ j ];

		int tval =  radialPlanePairs[ i ][ 0 ];
		radialPlanePairs[ i ][ 0 ] = radialPlanePairs[ minInd ][ 0 ];
		radialPlanePairs[ minInd ][ 0 ] = tval;
		tval = radialPlanePairs[ i ][ 1 ];
		radialPlanePairs[ i ][ 1 ] = radialPlanePairs[ minInd ][ 1 ];
		radialPlanePairs[ minInd ][ 1 ] = tval;	
	}

	//////////////////////////////////////////////////////////////////////////
	//debug, check if they are ordered in increasing orders
	/*for( int i = 0; i < size; i++ )
	{
		for( int j = 0; j < interPtIndVec[ i ].size(); j ++ )
		{
			cout<<interPtIndVec[ i ][ j ]<<" ";
		}
		cout<<"planes"<<radialPlanePairs[ i ][ 0 ]<<","<<radialPlanePairs[ i ][ 1 ]<<endl;		
	}*/
	//////////////////////////////////////////////////////////////////////////

	//step3, based on the sorting, 
	//merge those radial set of planes based on the intersection point indices
	int start, end;
	start = 0;	
	while( start < size )
	{
		end = start + 1;

		//search until the vec at position end is different from the one at start
		while( end < size && 
			compare2Vector( interPtIndVec[ start ], interPtIndVec[ end ] ) == 0	)
			end ++;

		//merge only there exists intersection points!
		if( interPtIndVec[ start ].size() <= 1 )	//at least 2 interpts.
		{
			start = end;
			continue;
		}

		//merge all the planes from pairs at positions in the range [start, end) 
		set<int> pSet;
		for( int i = start; i < end; i ++ )
		{
			pSet.insert( radialPlanePairs[ i ][ 0 ]  );
			pSet.insert( radialPlanePairs[ i ][ 1 ]  );
		}

		//only radial set with size >= 3 is considered, otherwise, general case, no need to 
		//handle differently.
		if( pSet.size() >= 3)
		{
			_gloabl_radialPlaneSet.push_back( pSet );
		}
		
		start = end;
	}	

	//clear the pairs, useless
	for( int i= 0; i< size; i ++) 
		radialPlanePairs[ i ].clear();
	radialPlanePairs.clear();
	for (unsigned int i = 0; i < interPtIndVec.size(); i++)
		interPtIndVec[ i ].clear();
	interPtIndVec.clear();	

	//////////////////////////////////////////////////////////////////////////
	//debug information, output all the sets of planes
	/*for( int i = 0; i < _gloabl_radialPlaneSet.size(); i++ )
	{
		set<int>::iterator iter = _gloabl_radialPlaneSet[ i ].begin();
		while( iter != _gloabl_radialPlaneSet[ i ].end() )
		{
			cout << *iter<<" ";
			iter ++;
		}
		cout<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
}

void resavePlanes( floatvector& param,	
				  vector<floatvector>& ctrvers,vector<intvector>& ctredges,					
				  vector<intvector>& ver2planelistpos, vector<intvector>& ver2planelist)
{
	int planeNum = graph_plane_edges.size();

	//step 0. plane parameters
	param.resize( 4 * planeNum );
	memcpy( &param[ 0 ], &graph_plane_params[ 0 ], sizeof( float ) * 4 * planeNum );

	ctrvers.resize( planeNum );
	ctredges.resize( planeNum );
	ver2planelistpos.resize( planeNum );
	for(int i = 0; i < planeNum; i ++ )
	{
		//step 1. go through each plane, save the index 
		int verNum = NodeIndex_Plane2Graph[ i ].size();
		ctrvers[ i ].resize( verNum * 3 );
		int pos = 0;
		for( int j = 0; j < verNum; j ++ )
		{
			int nodeInd = NodeIndex_Plane2Graph[ i ][ j ];
			memcpy( &ctrvers[ i ][ pos ], &graph_nodes[ 3 * nodeInd ] , sizeof( float ) * 3 );			
			pos += 3;
		}

		//step 2. go through the read in plane contour edges, resave them in the planes
		int edgeNum = graph_plane_edges[ i ].size() / 4;
		ctredges[ i ].resize( edgeNum * 4 );
		pos = 0;
		for( int j = 0; j < edgeNum; j++ )
		{
			for( int k = 0; k < 2; k ++ )	//set vertex index
			{
				ctredges[ i ][ pos ] = NodeIndex_Graph2Plane[ graph_plane_edges[ i ][ pos ] ][ i ];
				pos ++;
			}
			for( int k= 0; k < 2; k ++ )	//set material
			{
				ctredges[ i ][ pos ] =  graph_plane_edges[ i ][ pos ];
				pos ++;
			}
		}

		//step 3. go through each plane vertex, find the crspding graph nodes ,find out all the 
		//planes it is on
		ver2planelistpos[ i ].resize( verNum, -1 );
		for( int j = 0; j < verNum ;j ++ )
		{
			//check the vertex  j on plane i, find crspding node in the graph
			//from graph nodes to the existence map, find all the planes it's on
			int nodeI = NodeIndex_Plane2Graph[ i ][ j ];
			vector<int> planes;
			for( int k = 0;  k < planeNum; k ++ )
			{
				if( k == i )	//the current plane doesn't count
					continue;
				if( ExistGraphNodeInPlane[ nodeI ][ k ] )
				{
					planes.push_back( k );
				}
			}
			if( planes.size() == 0 )	//not an intersection point
				continue;

			//push the intersection plane list into the ver2planelist, and set pos for it
			ver2planelistpos[ i ][ j ] = ver2planelist.size();
			ver2planelist.push_back( planes );			
		}
	}

	//////////////////////////////////////////////////////////////////////////
	/*cout<<"ver2planelist:"<<endl;
	for( int  i = 0; i < planeNum; i++ )
	{
		cout<<"plane "<<i<<endl;
		int verNum = ver2planelistpos[ i ].size();
		for( int j= 0; j < verNum; j++ )
		{
			int pos = ver2planelistpos[ i ][ j ];
			if( pos == -1 )
				continue;

			cout<<"node in graph:"<<NodeIndex_Plane2Graph[ i ][ j ]<<":";
			for( int k= 0; k < ver2planelist[ pos ].size(); k ++ )
			{				
				cout<<ver2planelist[ pos ][ k] <<" ";
			}
			cout<<endl;
		}
	}*/
	//////////////////////////////////////////////////////////////////////////
}

void processCtrGraph( floatvector& param,	
					 vector<floatvector>& ctrvers,vector<intvector>& ctredges,					
					 vector<intvector>& ver2planelistpos, vector<intvector>& ver2planelist)
{
	//step1. mark existence
	markExist();

	//step2. get the radial sets of planes that share one common line
	radialSetAnalyze();

	//step3. resave them into the conventional format, set the ver2planelistpos and ver2planelist
	resavePlanes(param, ctrvers, ctredges, ver2planelistpos, ver2planelist );

	//step 4. clean all the temp variables in this file
	graph_nodes.clear();
	for (unsigned int i = 0; i < graph_plane_edges.size(); i++)
	{
		graph_plane_edges[ i ].clear();		
	}
	graph_plane_params.clear();
	graph_plane_edges.clear();

	for (unsigned int i = 0; i < ExistGraphNodeInPlane.size(); i++)
		ExistGraphNodeInPlane[ i ].clear();
	ExistGraphNodeInPlane.clear();

	for (unsigned int i = 0; i < NodeIndex_Graph2Plane.size(); i++)
		NodeIndex_Graph2Plane[ i ].clear();
	NodeIndex_Graph2Plane.clear();
    
	for(unsigned int i = 0; i < NodeIndex_Plane2Graph.size(); i++ )
		NodeIndex_Plane2Graph[  i ].clear();
	NodeIndex_Plane2Graph.clear();
}

void ContourHandler::readContour3_v3(
									 const int filenum, char** filenames, 
									 floatvector& param,	
									 vector<floatvector>& ctrvers,vector<intvector>& ctredges,
									 float bbox[ 6 ], bool& bboxset,
									 vector<intvector>& ver2planelistpos,
									 vector<intvector>& ver2planelist	 )
{
	cout<<"Start reading "<<filenames[ 0 ]<<"..... ";

	//later on, i decide that multiple files loading is not allowable.
	//only one file at a time. so when coming here, one file is guaranteed.
	//only call the following funcion once, insteand of filenum times.
	//read file one by one
	readOneContour3_v3(filenames[ 0 ],bbox,	bboxset);

	//check if success
    if( graph_nodes.size() == 0 )	//not successful!
	{		
		return;
	}

	//cout<<"Reading contour3 successfully!"<<endl;
	
	//process the input graph to make them the save as old formats.
	processCtrGraph(param, ctrvers, ctredges, ver2planelistpos, ver2planelist);

	//normalize the normals of the planes
	int planenum = param.size()/4;
	for( int i = 0; i < planenum; i ++)
	{
		float veclen = MyMath::vectorLen( param[ i * 4 ], param[ i*4 + 1 ], param[ i*4+ 2 ]);
		for(int j = 0; j < 4; j ++)
			param[ 4*i + j] = param[ 4*i + j]/veclen;
	}
	cout<<"Done!\n";
}
