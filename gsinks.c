/* gsinks.c */

/* 
Reads input files of digraphs and counts the number of unlabeled digraphs with a single global sink that can be created by adding a new vertex and adding edges from the digraph to the new vertex. 

Based heavily on the following:
nauty version 2.7 
especially vcolg.c 

Also, the tarjan algorithm for partitioning digraph into strongly-connected components (SCCs)

Input data from B. D McKay:
https://users.cecs.anu.edu.au/~bdm/data/digraphs.html as of 2021-07-02


Basic logic:

for each graph in the input file

Run the tarjan algorithm to partition the digraph into strongly-connected components (SCCs). Along the way, identify which SCCs are leaves of the resulting DAG of SCCs.  

Then, use code from vcolg.c to identify two-colorings of the digraph.

Filter out any colorings that don't color at least one vertex of each Leaf SCC.

Count the remaining colorings. If generating graphs, create a new vertex (which will become the global sink), and extend new edges to that new vertex from every vertex colored '1'. 


Building instructions
SEE readme.MD

*/

#define USAGE \
  "gsinks [opts] N"

#define HELPTEXT \
" gsinks : generate and count digraphs with one global sink\n\
\n\
     -q     don't show total count for the number of vertices\n\
     -d     generate, and show in d6 formmat \n\
     -l     self-loops allowed \n\
     N      vertex count  (default: start at 1 and go up)\n\
"

/* Nauty-required definitions before any includes */
#define MAXN 32 
#define WORDSIZE 32

#include "gtools.h"
#include "naugroup.h"
#include "nautinv.h"


// From vcolg.c
static int col[MAXN];
static boolean first;
static int lastreject[MAXN];
static boolean lastrejok;
static unsigned long groupsize;
static unsigned long newgroupsize;
static int fail_level;

// switches
static boolean qswitch, dswitch, lswitch;

// for counting and generating digraphs with one global sink
static long long totalCount = 0;

// for tarjan algorithm 
// trajan breaks graphs into strongly connected components (SCCs)
// and the connections between SCCs is a directed acyclic graph (DAG)
struct sccinfo {setword sccVertices; boolean isLeaf; int sccSize;};
struct sccinfo sccinfos[MAXN];
int currentSCC;
boolean isFirstSCC;
struct vinfo {int index; int lowlink; boolean onstack; set descendents;};
struct vinfo vinfos[MAXN];
static int vertex_index = 0;
#define UNDEFINED (-1)

// a very simple integer stack implementation, used by tarjan
static int stack[MAXN];     
static int stack_top = -1;  
static int pop() {
    int data = stack[stack_top];
    stack_top--;   
    return data;
}
static int push(int data) {
    stack_top++;   
    stack[stack_top] = data;
}



// A callout from colour_digraph code, as digraphs get colored.  
// Every node colored '1' will be connected to the new global sink vertex.
// To ensure a single gloabl sink, there must be at least one connection 
// from every leaf SCC to the new global sink vertex. This function
// rejects the graph colorings that don't fulfill  that requirement, and, 
// if requested, contructs and outputs every graph that meets the requirement. 

void filter_and_output(graph*,int*,int,int);
void filter_and_output(graph* g,int* col,int m,int n){
    set * gi;
    int i,j;
    boolean leafcoloured;
    
    // reject graphs that don't have at least one '1' in every SCC leaf
    for (i=0; i<currentSCC; i++)
    {
        if (!sccinfos[i].isLeaf) continue;
        
        leafcoloured = FALSE;
        j = nextelement(&sccinfos[i].sccVertices,m,UNDEFINED);
        do
        {
            if (col[j]) {
                leafcoloured = TRUE;
                break;
            }
            j = nextelement(&sccinfos[i].sccVertices,m,j);
        } while (j != UNDEFINED);
        if (FALSE == leafcoloured) return;
    }
   
    // Now create the actual single-sink digraphs by adding a new edge & connecting all colored edges to it. 
    if (dswitch)
    {
        graph gnew [MAXN*MAXM];
        EMPTYGRAPH(gnew,m,n+1);
        memcpy(gnew, g, n * sizeof(set));
        
        for (j=0; j<n; j++)
        {
            if (col[j]){
                gi = GRAPHROW(gnew,j,m);
                ADDELEMENT(gi,n);    
            }
        }
        writed6(stdout, gnew, m, n+1);  
    }
    totalCount++;
}

/**************************************************************************/
//Code from vcolg.c, simplified where possible


static int
ismax(int *p, int n)
/* test if col^p <= col */
{
    int i,k;
    int fail;

    fail = 0;
    for (i = 0; i < n; ++i)
    {
	k = p[i];
	if (k > fail) fail = k;
        if (col[k] > col[i])
        {
            fail_level = fail;
            return FALSE;
        }
        else if (col[k] < col[i]) return TRUE;
    }

    ++newgroupsize;
    return TRUE;
}

/**************************************************************************/

static void
testmax(int *p, int n, int *abort)
/* Called by allgroup2. */
{
    int i;

    if (first)
    {                       /* only the identity */
        first = FALSE;
        return;
    }

    if (!ismax(p,n))
    {
        *abort = 1;
        for (i = 0; i < n; ++i) lastreject[i] = p[i];
        lastrejok = TRUE;
    }
}

/**************************************************************************/

static int
trythisone(grouprec *group, graph *g, boolean digraph, int m, int n)
/* Try one solution, accept if maximal. */
/* Return value is level to return to. */
{
    int i,j;
    boolean accept;
    graph *gi;
    size_t ne;

    newgroupsize = 1;

    if (!group || groupsize == 1)
        accept = TRUE;
    else if (lastrejok && !ismax(lastreject,n))
        accept = FALSE;
    else if (lastrejok && groupsize == 2)
        accept = TRUE;
    else
    {
        newgroupsize = 1;
        first = TRUE;

        if (allgroup2(group,testmax) == 0)
            accept = TRUE;
        else
            accept = FALSE;
    }

    if (accept)
    {        
        filter_and_output(g,col,m,n);

        return n-1;
    }
    else
        return fail_level-1;
}

/**************************************************************************/

static int
scan(int level, graph *g, boolean digraph, int *prev, long minedges, long maxedges,
    long sofar, long numcols, grouprec *group, int m, int n)
/* Recursive scan for default case */
/* Returned value is level to return to. */
{
    int left;
    long min,max,k,ret;

    if (level == n)
        return trythisone(group,g,digraph,m,n);

    left = n - level - 1;
    min = minedges - sofar - numcols*left;
    if (min < 0) min = 0;
    max = maxedges - sofar;
    if (max >= numcols) max = numcols - 1;
    if (prev[level] >= 0 && col[prev[level]] < max)
	max = col[prev[level]];

    for (k = min; k <= max; ++k)
    {
        col[level] = k;
        ret = scan(level+1,g,digraph,prev,minedges,maxedges,sofar+k,numcols,group,m,n);
	if (ret < level) return ret;
    }

    return level-1;
}

/**************************************************************************/

static void
colourdigraph(graph *g, int nfixed, long minedges, long maxedges,
         long numcols, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[MAXN];
    grouprec *group;
    int i,j,k,nloops;
    size_t ii;
    set *gi,*gj,*gci,*gcj;
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    boolean loop[MAXN];
    int prev[MAXN]; /* If >= 0, earlier point that must have greater colour */
    int weight[MAXN];
    int region,start,stop;
    DYNALLSTAT(graph,gconv,gconv_sz);

    if (n > MAXN) gt_abort(">E vcolg: MAXN exceeded\n");
    

    DYNALLOC2(graph,gconv,gconv_sz,n,m,"colourdigraph");

    nloops = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i))
	{
	    DELELEMENT(gi,i);
	    loop[i] = TRUE;
	    ++nloops;
	}
	else
	    loop[i] = FALSE;

    for (ii = 0; ii < m*(size_t)n; ++ii) gconv[ii] = g[ii];
    converse(gconv,m,n);

    for (region = 0; region < 2; ++region)
    {
	if (region == 0)
	{
	    if (nfixed == 0) continue;
	    start = 0;
	    stop = nfixed;
	    if (stop > n) stop = n;
	}
	else
	{
	    if (nfixed >= n) continue;
	    start = nfixed;
	    stop = n;
	}
	
	for (i = start,
                    gi = g + m*(size_t)start, gci = gconv + m*(size_t)start;
             i < stop; ++i, gi += m, gci += m)
	{
            /* Find most recent equivalent j. */
	    for (j = i-1, gj = gi-m, gcj = gci-m; j >= start;
                                                   --j, gj -= m, gcj -= m)
	    {
		if (loop[j] != loop[i]
                       || ISELEMENT(gi,j) != ISELEMENT(gj,i)) continue;
		for (k = 0; k < m; ++k)
                     if (gi[k] != gj[k] || gci[k] != gcj[k]) break;
		if (k < m)
		{
		    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
		    FLIPELEMENT(gci,i); FLIPELEMENT(gcj,j);
		    for (k = 0; k < m; ++k)
                        if (gi[k] != gj[k] || gci[k] != gcj[k]) break;
		    FLIPELEMENT(gci,i); FLIPELEMENT(gcj,j);
		    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
		}
		if (k == m) break;
	    }
	    if (j >= start)
	    {
		prev[i] = j;
		weight[i] = weight[j] + 1;
	    }
	    else
	    {
		prev[i] = -1;
		weight[i] = 0;
	    }
	}
    }

    for (i = nfixed; i < n; ++i) weight[i] += nfixed;

    if (maxedges == NOLIMIT || maxedges > n*numcols) maxedges = n*numcols;
    if (n*numcols < minedges) return;

    if (n == 0)
    {
        scan(0,g,TRUE,prev,minedges,maxedges,0,numcols,FALSE,m,n);
        return;
    }

    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    options.defaultptn = FALSE;
    options.digraph = TRUE;
    options.invarproc = adjacencies;
    options.maxinvarlevel = n;

    setlabptn(weight,lab,ptn,n);

    if (nloops > 0)
        for (i = 0, gi = g; i < n; ++i, gi += m)
	    if (loop[i]) ADDELEMENT(gi,i);
 
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,MAXN,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    if (stats.numorbits < n)
    {
	j = n;
	for (i = 0; i < n; ++i)
	    if (orbits[i] < i && orbits[i] < j) j = orbits[i];

	for (i = j + 1; i < n; ++i)
	    if (orbits[i] == j) prev[i] = j;
    }

    lastrejok = FALSE;
    for (i = 0; i < n; ++i) col[i] = 0;

    scan(0,g,TRUE,prev,minedges,maxedges,0,numcols,group,m,n);
}

/***************************/
// TARJAN algorithm
/*********************/

void tarjan(graph * g, int m, int n);
void strongconnect(graph * g,int v,int n,int m);

void tarjan(graph * g, int m, int n)
{
    int i,j;
     
    stack_top = -1;
    currentSCC = 0;
    isFirstSCC = TRUE;
    
    for (i=0; i<n; i++){ 
        vinfos[i].index = UNDEFINED;
        vinfos[i].lowlink = UNDEFINED;
        vinfos[i].onstack = FALSE;
        EMPTYSET(&vinfos[i].descendents,m);
        
        EMPTYSET(&sccinfos[i].sccVertices,m);
        sccinfos[i].isLeaf = FALSE;
        sccinfos[i].sccSize = 0;
    }
    
    
    for (int v=0; v<n; v++)
    {
        if (UNDEFINED == vinfos[v].index){
            strongconnect(g,v,n,m);
        }
    }
    for (i=0; i<currentSCC; i++)
    {
        j = nextelement(&sccinfos[i].sccVertices,m,UNDEFINED);
        do
        {
            j = nextelement(&sccinfos[i].sccVertices,m,j);
        } while (j != UNDEFINED);
        if (sccinfos[i].isLeaf)
        {
        }
    }    
    
}

void strongconnect(graph * g,int v,int n,int m)
{
    set * gp;
    int w = UNDEFINED;
    setword descendents;
    int sccSize;
    
    vinfos[v].index = vertex_index;
    vinfos[v].lowlink = vertex_index;
    vertex_index++;
    push(v);
    vinfos[v].onstack = TRUE;
    
    gp = GRAPHROW(g,v,m);
    w = nextelement(gp,m,w);
    
    for (; w > UNDEFINED; w = nextelement(gp,m,w))
    {
        ADDELEMENT(&vinfos[v].descendents,w);
        if (vinfos[w].index == UNDEFINED)
        {
            strongconnect(g,w,n,m);
            
            if (vinfos[v].lowlink > vinfos[w].lowlink)
            {
                vinfos[v].lowlink  =  vinfos[w].lowlink;
            }
        }
        else if (vinfos[w].onstack)
        {
            if (vinfos[v].lowlink > vinfos[w].index)
            {
                vinfos[v].lowlink = vinfos[w].index;
            }
        }
    }  
    EMPTYSET(&descendents, m);
    
    if (vinfos[v].lowlink == vinfos[v].index){
        sccSize=0;
        // start a new strongly connected component
        do {
            w = pop();
            vinfos[w].onstack = FALSE;
            sccSize++;
            // add w to current strongly connected component
            ADDELEMENT(&sccinfos[currentSCC].sccVertices, w);
            // accumulate descendents to check for leafiness
            ADDELEMENT(&descendents, w);
            descendents = UNION(descendents,vinfos[w].descendents);
        } while (w != v);
        sccinfos[currentSCC].sccSize = sccSize;
        
        
        if (isFirstSCC){
            sccinfos[currentSCC].isLeaf = TRUE;
            isFirstSCC = FALSE;
        }
        else if (0 == SETDIFF(descendents,sccinfos[currentSCC].sccVertices)) {
            sccinfos[currentSCC].isLeaf = TRUE;
        }      
            
        // done with recording this strongly connected component
        currentSCC++;
    }
}
            
/**********************************************************************/
/* Top-level function: read input arguments, open appropriate 
/* input files & launch the algorithm */
/**********************************************************************/

#define INFILE_PREFIX "dig"
#define INFILE_LOOP_MODIFIER 'l'
#define INFILE_SUFFIX ".d6"

int
main(int argc, char *argv[])
{
    graph *g;
    int m,n,codetype;
    char infilename[10];    // dig[l][n].d6
    FILE *infile;
    int i,j;
    char *arg;
    char * filepart;
    char * endptr;
    boolean badargs,digraph;
    long countN = 0;
    int startN, endN;

    HELP; PUTVERSION;
    
    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    // default values
    qswitch = FALSE;
    dswitch = FALSE;
    lswitch = FALSE;

    infilename[0] = '\0';
    
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-')
        {
            switch (arg[1]) {
                case 'q': qswitch = TRUE; break;
                case 'd': dswitch = TRUE; break;
                case 'l': lswitch = TRUE; break;
                default: badargs = TRUE;    
            }
        }
        else
        {
            countN = strtol(arg,&endptr,10);
            if (countN <= 0){badargs = TRUE;}
        }
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (countN == 1) 
    {
        // Need to write special code for this case 
    }
    if (countN)
    {
        startN = countN-1;
        endN = countN;
    }
    else
    {
        startN = 1;
        endN = 10;  // Arbitraily set at 10. It will curently stop at less because of data file limits; 
    }
    
    
    for (i = startN; i < endN; i++)
    {
        // contruct the input filename to be used
        strcpy(infilename,INFILE_PREFIX);
        filepart = infilename + strlen(INFILE_PREFIX);
        if (lswitch){
            *filepart++ = INFILE_LOOP_MODIFIER;
        }
        snprintf(filepart, filepart - infilename, "%d", i);
        strcat(infilename, INFILE_SUFFIX);
        
        // process each graph in each input file
        infile = opengraphfile(infilename,&codetype,FALSE,1);
        if (!infile) exit(1);
        
        while (NULL != (g = readgg_inc(infile,NULL,0,&m,&n,
                            NULL,1,1,&digraph)))
        {
            tarjan(g, m, n);
            // Now color each graph, now that we know the SCCs 
            colourdigraph(g,0,0,NOLIMIT,2,m,n);
            // colourdigraph will call out to filter and count the graphs and output them if requested
        }
        if (!qswitch){
            fprintf(stderr,"%lld\n",totalCount);
        }
        totalCount = 0;
    }
    exit(0);
}
