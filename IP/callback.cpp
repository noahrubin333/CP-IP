/* Custom callback function for Gurobi that generates clique cuts in the MOLS problem */

// Generate cuts from cliques of type II
#define CLIQUETYPEII
// Generate cuts from cliques of type III
#define CLIQUETYPEIII

int totalcuts = 0;

#include <tuple>
#include <vector>
typedef tuple<int, int, int, int> quadruple;

// Return a vector of quadruples that have exactly two entries in common with the quadruple (i,j,k,l)
vector<quadruple> intersection_2(int i, int j, int k, int l)
{
    vector<quadruple> X = {};

    for(int ip=0; ip<N; ip++)
        for(int jp=0; jp<N; jp++)
            if(ip != i && jp != j)
                X.push_back(make_tuple(ip,jp,k,l));
    for(int ip=0; ip<N; ip++)
        for(int kp=0; kp<N; kp++)
            if(ip != i && kp != k)
                X.push_back(make_tuple(ip,j,kp,l));
    for(int ip=0; ip<N; ip++)
        for(int lp=0; lp<N; lp++)
            if(ip != i && lp != k)
                X.push_back(make_tuple(ip,j,k,lp));
    for(int jp=0; jp<N; jp++)
        for(int kp=0; kp<N; kp++)
            if(jp != i && kp != k)
                X.push_back(make_tuple(i,jp,kp,l));
    for(int jp=0; jp<N; jp++)
        for(int lp=0; lp<N; lp++)
            if(jp != i && lp != k)
                X.push_back(make_tuple(i,jp,k,lp));
    for(int kp=0; kp<N; kp++)
        for(int lp=0; lp<N; lp++)
            if(kp != i && lp != k)
                X.push_back(make_tuple(i,j,kp,lp));

    return X;
}

// Return a vector of quadruples that have exactly three entries in common with the quadruple (i,j,k,l)
vector<quadruple> intersection_3(int i, int j, int k, int l)
{
    vector<quadruple> X = {};

    for(int ip=0; ip<N; ip++)
        if(ip != i)
            X.push_back(make_tuple(ip,j,k,l));
    for(int jp=0; jp<N; jp++)
        if(jp != j)
            X.push_back(make_tuple(i,jp,k,l));
    for(int kp=0; kp<N; kp++)
        if(kp != k)
            X.push_back(make_tuple(i,j,kp,l));
    for(int lp=0; lp<N; lp++)
        if(lp != l)
            X.push_back(make_tuple(i,j,k,lp));

    return X;
}

// Macro that returns true when exactly one of the arguments 'x' is equal to the argument 'xp'
#define INT1(i,j,k,l,ip,jp,kp,lp) ((i==ip && j!=jp && k!=kp && l!=lp)||(i!=ip && j==jp && k!=kp && l!=lp)||(i!=ip && j!=jp && k==kp && l!=lp)||(i!=ip && j!=jp && k!=kp && l==lp))

#ifdef VERBOSE
    #define LOGCUTII cout << "Added custom cut derived from clique of type II" << endl
    #define LOGCUTIII cout << "Added custom cut derived from clique of type III" << endl
#else
    #define LOGCUTII
    #define LOGCUTIII
#endif

class MyCallback: public GRBCallback
{
    void callback()
    {
#ifdef VERBOSE
        if(where==GRB_CB_MIP)
        {
            double nodenum = getDoubleInfo(GRB_CB_MIP_NODCNT);
            int cuts = getIntInfo(GRB_CB_MIP_CUTCNT);
            cout << "Node " << nodenum << " Applied Cuts " << cuts << " Generated Cuts " << totalcuts << endl;
        }
#endif
        if(where==GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
        {
#ifdef CLIQUETYPEII
            // Check for and add cuts generated from cliques of type II
            const int v = 5;
            double d[N][N][N][N] = {};
            for(int i=0; i<N; i++)
            {   for(int j=0; j<N; j++)
                {   for(int k=0; k<N; k++)
                    {   for(int l=0; l<N; l++)
                        {   if(getNodeRel(xvar[i][j][k][l]) >= 1/(double)(v*N))
                                d[i][j][k][l] += getNodeRel(xvar[i][j][k][l]);
                            if(d[i][j][k][l] > 1)
                            {
                                GRBLinExpr expr = 0;
                                for(int ip=0; ip<N; ip++)
                                {
                                    if(ip != i)
                                        expr += xvar[ip][j][k][l];
                                }
                                for(int jp=0; jp<N; jp++)
                                {
                                    if(jp != j)
                                        expr += xvar[i][jp][k][l];
                                }
                                for(int kp=0; kp<N; kp++)
                                {
                                    if(kp != k)
                                        expr += xvar[i][j][kp][l];
                                }
                                for(int lp=0; lp<N; lp++)
                                {
                                    if(lp != l)
                                        expr += xvar[i][j][k][lp];
                                }
                                addCut(expr, GRB_LESS_EQUAL, 1.0);
                                totalcuts++;
                                LOGCUTII;
                            }
                        }
                    }
                }
            }
            for(int i=0; i<N; i++)
            {   for(int j=0; j<N; j++)
                {   for(int k=0; k<N; k++)
                    {   for(int l=0; l<N; l++)
                        {   if(d[i][j][k][l] > (v-4)/(double)v)
                            {
                                double sum = 0;
                                for(int ip=0; ip<N; ip++)
                                {
                                    if(ip != i)
                                        sum += getNodeRel(xvar[ip][j][k][l]);
                                }
                                for(int jp=0; jp<N; jp++)
                                {
                                    if(jp != j)
                                        sum += getNodeRel(xvar[i][jp][k][l]);
                                }
                                for(int kp=0; kp<N; kp++)
                                {
                                    if(kp != k)
                                        sum += getNodeRel(xvar[i][j][kp][l]);
                                }
                                for(int lp=0; lp<N; lp++)
                                {
                                    if(lp != l)
                                        sum += getNodeRel(xvar[i][j][k][lp]);
                                }

                                if(sum > 1)
                                {
                                    GRBLinExpr expr = 0;
                                    for(int ip=0; ip<N; ip++)
                                    {
                                        if(ip != i)
                                            expr += xvar[ip][j][k][l];
                                    }
                                    for(int jp=0; jp<N; jp++)
                                    {
                                        if(jp != j)
                                            expr += xvar[i][jp][k][l];
                                    }
                                    for(int kp=0; kp<N; kp++)
                                    {
                                        if(kp != k)
                                            expr += xvar[i][j][kp][l];
                                    }
                                    for(int lp=0; lp<N; lp++)
                                    {
                                        if(lp != l)
                                            expr += xvar[i][j][k][lp];
                                    }
                                    addCut(expr, GRB_LESS_EQUAL, 1.0);
                                    totalcuts++;
                                    LOGCUTII;
                                }
                            }
                        }
                    }
                }
            }
#endif
#ifdef CLIQUETYPEIII
            // Check for and add cuts generated from cliques of type III
            for(int i=0; i<N; i++)
            {   for(int j=0; j<N; j++)
                {   for(int k=0; k<N; k++)
                    {   for(int l=0; l<N; l++)
                        {
                            double xc = getNodeRel(xvar[i][j][k][l]);
                            if(xc > 0.25 && xc < 1)
                            {
                                vector<quadruple> X = intersection_2(i,j,k,l);
                                int xsize = X.size();
                                for(int s=0; s<xsize; s++)
                                {
                                    int ip = get<0>(X[s]);
                                    int jp = get<1>(X[s]);
                                    int kp = get<2>(X[s]);
                                    int lp = get<3>(X[s]);
                                    double xt = getNodeRel(xvar[ip][jp][kp][lp]);
                                    if(xt > (1-xc)/(double)3)
                                    {
                                        vector<quadruple> Y = intersection_3(ip,jp,kp,lp);
                                        int ysize = Y.size();
                                        for(int s2=0; s2<ysize; s2++)
                                        {
                                            int ipp = get<0>(Y[s2]);
                                            int jpp = get<1>(Y[s2]);
                                            int kpp = get<2>(Y[s2]);
                                            int lpp = get<3>(Y[s2]);
                                            if(INT1(i,j,k,l,ipp,jpp,kpp,lpp))
                                            {
                                                if(i==ipp && xc+getNodeRel(xvar[i][j][kp][lp])+getNodeRel(xvar[i][jp][k][lp])+getNodeRel(xvar[i][jp][kp][l])>1)
                                                {   LOGCUTIII; totalcuts++; addCut(xvar[i][j][k][l]+xvar[i][j][kp][lp]+xvar[i][jp][k][lp]+xvar[i][jp][kp][l], GRB_LESS_EQUAL, 1.0); }
                                                if(j==jpp && xc+getNodeRel(xvar[i][j][kp][lp])+getNodeRel(xvar[ip][j][k][lp])+getNodeRel(xvar[ip][j][kp][l])>1)
                                                {   LOGCUTIII; totalcuts++; addCut(xvar[i][j][k][l]+xvar[i][j][kp][lp]+xvar[ip][j][k][lp]+xvar[ip][j][kp][l], GRB_LESS_EQUAL, 1.0); }
                                                if(k==kpp && xc+getNodeRel(xvar[i][jp][k][lp])+getNodeRel(xvar[ip][j][k][lp])+getNodeRel(xvar[ip][jp][k][l])>1)
                                                {   LOGCUTIII; totalcuts++; addCut(xvar[i][j][k][l]+xvar[i][jp][k][lp]+xvar[ip][j][k][lp]+xvar[ip][jp][k][l], GRB_LESS_EQUAL, 1.0); }
                                                if(l==lpp && xc+getNodeRel(xvar[i][jp][kp][l])+getNodeRel(xvar[ip][j][kp][l])+getNodeRel(xvar[ip][jp][k][l])>1)
                                                {   LOGCUTIII; totalcuts++; addCut(xvar[i][j][k][l]+xvar[i][jp][kp][l]+xvar[ip][j][kp][l]+xvar[ip][jp][k][l], GRB_LESS_EQUAL, 1.0); }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
#endif
        }
    }
};

