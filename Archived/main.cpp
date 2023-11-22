#include <thread>
#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <papi.h>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>
#include <stack>
#include <queue>
#include <utility>
#include <mutex>
#include <string>
#include "metis.h"

unsigned int INF = -1;
// std::atomic<unsigned int> globalThreadCall[3];

void pathDouble(unsigned int nVertex, unsigned int* apsp);

constexpr unsigned char SUPERNODE_ND_BFS = 0b00000001;
constexpr unsigned char SUPERNODE_ND_MET = 0b00000010;

class Supernode{//nested dissection result
public:
    Supernode(unsigned char supernodeFlag):supernodeFlag(supernodeFlag){
        parent = nullptr;
        VertexIdx = nullptr;
        length = 0;
        C1 = nullptr;
        C2 = nullptr;
    }
    ~Supernode(){
        if(supernodeFlag & SUPERNODE_ND_BFS){
            if(parent == nullptr)
                delete[] VertexIdx;
            if(C1 != nullptr)
                delete C1;
            if(C2 != nullptr)
                delete C2;
        }
        else if(supernodeFlag & SUPERNODE_ND_MET){
            if(parent == nullptr){
                freeVertex();
            }
            if(C1 != nullptr)
                delete C1;
            if(C2 != nullptr)
                delete C2;
        }
    }
    void freeVertex(){
        if(C1 == nullptr){
            delete[] VertexIdx;
        }
        else{
            C1->freeVertex();
        }
    }
    unsigned char supernodeFlag;
    unsigned int* VertexIdx;
    unsigned int length;
    Supernode* C1; //If nested use this
    Supernode* C2; //If nested use this
    Supernode* parent;
    // It forms E-tree naturally
};

void printBFS(Supernode* root){
    std::queue<std::pair<Supernode *, int>> bfsqueue;
    bfsqueue.emplace(root, 0);
    while(!bfsqueue.empty()){
        printf("[%d(%d)] %p %p %p %p)",bfsqueue.front().second,bfsqueue.front().first->length,bfsqueue.front().first,bfsqueue.front().first->parent,bfsqueue.front().first->C1,bfsqueue.front().first->C2);
        unsigned int from = 0;
        for(int i = 0; i < bfsqueue.front().first->length; ++i){
            printf(" %d",bfsqueue.front().first->VertexIdx[i]);
        }printf("\n");
        if(bfsqueue.front().first->C1 != nullptr)
            bfsqueue.emplace(bfsqueue.front().first->C1, bfsqueue.front().second + 1);
        if(bfsqueue.front().first->C2 != nullptr)
            bfsqueue.emplace(bfsqueue.front().first->C2, bfsqueue.front().second + 1);
        bfsqueue.pop();
    }

}

unsigned int countDepth(Supernode* node){
    if(node == nullptr)
        return 0;
    return std::max(countDepth(node->C1) + 1, countDepth(node->C2) + 1);
}

unsigned int countNum(Supernode* node){
    if(node == nullptr)
        return 0;
    return countNum(node->C1) + countNum(node->C2) + 1;
}

unsigned int countVertexNumUp(Supernode* node){
    if(node == nullptr)
        return 0;
    return countVertexNumUp(node->parent) + node->length;
}

unsigned int countVertexNumDown(Supernode* node){
    if(node == nullptr)
        return 0;
    return countVertexNumDown(node->C1) + 
            countVertexNumDown(node->C2) + node->length;
}

unsigned int countEleNum(Supernode* root){
    return countVertexNumDown(root);
}

void initResult(int nVertex, unsigned int** distance, unsigned int*** apsp){
    *apsp = new unsigned int*[nVertex];

    unsigned int** apsp_temp = *apsp;

    for(int i = 0; i < nVertex; ++i){
        apsp_temp[i] = new unsigned int[nVertex];
        for(int j = 0; j < nVertex; ++j){
            apsp_temp[i][j] = distance[i][j];
        }
    }
}

void freeResult(int nVertex, unsigned int** apsp){
    for(int i = 0; i < nVertex; ++i){
        delete[] apsp[i];
    }
    delete[] apsp;
}

unsigned int decomposeGPU(int nVertex, int i, int j = 0);

void initResultGPU(int nVertex, unsigned int** distance, unsigned int** apsp){
    *apsp = new unsigned int[nVertex * nVertex];

    unsigned int* apsp_temp = *apsp;

    for(int i = 0; i < nVertex; ++i){
        for(int j = 0; j < nVertex; ++j){
            apsp_temp[decomposeGPU(nVertex, i, j)] = distance[i][j];
        }
    }
}

void freeResultGPU(unsigned int* apsp){
    delete[] apsp;
}

void updateEdge(int i, int j, int k, int nVertex, unsigned int** apsp){
    if(apsp[i][k] == INF || apsp[k][j] == INF)
        return;
    if(apsp[i][j] > apsp[i][k] + apsp[k][j]){
        apsp[i][j] = apsp[i][k] + apsp[k][j];
    }
}

void updateEdge_atomic(int i, int j, int k, std::mutex** mutexBlock, int nVertex, unsigned int** apsp){
    if(apsp[i][k] == INF || apsp[k][j] == INF)
       return;
    if(apsp[i][j] > apsp[i][k] + apsp[k][j]){
        mutexBlock[i][j].lock();
        if(apsp[i][j] > apsp[i][k] + apsp[k][j]){
            apsp[i][j] = apsp[i][k] + apsp[k][j];
        }
        mutexBlock[i][j].unlock();  
    }
}

unsigned int getValue(unsigned int range, unsigned int size){
    // return the start value of range
    if(range == 0)
        return 0;
    if(range == 1)
        return size;
    int v_1 = int(range/2);
    int v_2 = int((range + 1)/2);
    return size*  (v_1 + v_2) - v_2 * (v_2 - 1) / 2 - v_1 * (v_1 + 1)/2;
}

unsigned int getRange(unsigned int value, unsigned int size){
    unsigned int low = 0;
    unsigned int high = 2 * size;
    unsigned int mid,midVS,midVE;
    while(true){
        mid = (low + high) / 2;
        midVS = getValue(mid, size);
        midVE = getValue(mid + 1, size);
        if(low > high)
            exit(0);
        if(midVS <= value && value < midVE){
            return mid;
        }
        if(midVS == value && value == midVE){
            return mid;
        }
        else if(value < midVS){
            high = mid - 1;
        }
        else{
            low = mid + 1;
        }
    }
}

void updateEdgePerElem(int i, int j, int k, int nVertex, unsigned int** apsp){
    updateEdge(i,j,k,nVertex,apsp);
}
void updateEdgePerLine(int j, int k, int nVertex, unsigned int** apsp){
    for(int i = 0; i < nVertex; ++i)
        updateEdge(i,j,k,nVertex,apsp);
}
void updateEdgePerThread(int threadIdx, int k, int nThread, int nVertex, unsigned int** apsp){
    int length = int(sqrt(nThread));
    if(threadIdx >= length * length)
        return;
    
    int unitWidth  = nVertex / length;
    int unitHeight = unitWidth;

    int threadIdxX = threadIdx % length;
    int threadIdxY = threadIdx / length;
    
    int fromi,toi,fromj,toj;
    fromi = threadIdxX * unitWidth;
    fromj = threadIdxY * unitHeight;
    
    if(threadIdxX == (length - 1))
        toi = nVertex;
    else
        toi = (threadIdxX + 1) * unitWidth;
        
    if(threadIdxY == (length - 1))
        toj = nVertex;
    else
        toj = (threadIdxY + 1) * unitHeight;

    for(int i = fromi; i < toi; ++i)
        for(int j = fromj; j < toj; ++j)
            updateEdge(i,j,k,nVertex,apsp);
}
void updateEdgeBlocked(unsigned int fromI, unsigned int toI, unsigned int fromJ, unsigned int toJ, unsigned int fromK, unsigned int toK, int nVertex, unsigned int** apsp){
    for(int k = fromK; k < toK; ++k)
        for(int i = fromI; i < toI; ++i)
            for(int j = fromJ; j < toJ; ++j)
                updateEdge(i,j,k,nVertex,apsp);
}

enum NodedMODE{
    COL,ROW,OUTTER
};

void updateEdgeNoded(unsigned int* vertexIdxI, unsigned int * vertexIdxJ, unsigned int * vertexIdxK,
            unsigned int fromI, unsigned int toI, unsigned int fromJ, unsigned int toJ, unsigned int fromK, unsigned int toK, 
            int nVertex, unsigned int** apsp){
    for(int k = fromK; k < toK; ++k)
            for(int i = fromI; i < toI; ++i)
                for(int j = fromJ; j < toJ; ++j)
                    updateEdge(vertexIdxI[i],vertexIdxJ[j],vertexIdxK[k],nVertex,apsp);
}

void SetUpVertexListDown(unsigned int* vertexIdx, Supernode* node, unsigned int& idxInUse){
    if(node == nullptr)
        return;
    for(int i = 0; i < node->length; ++i, ++idxInUse){
        vertexIdx[idxInUse] = node->VertexIdx[i];
    }
    SetUpVertexListDown(vertexIdx, node->C1, idxInUse);
    SetUpVertexListDown(vertexIdx, node->C2, idxInUse);
}

void SetUpVertexListUp(unsigned int* vertexIdx, Supernode* node, unsigned int& idxInUse){
    if(node == nullptr)
        return;
    for(int i = 0; i < node->length; ++i, ++idxInUse){
        vertexIdx[idxInUse] = node->VertexIdx[i];
    }
    SetUpVertexListUp(vertexIdx, node->parent, idxInUse);
}

void SetUpVertexList(unsigned int** vertexIdx, Supernode* node, int* nVertex, int* separatorFrom){
    int nAncester = countVertexNumUp(node->parent);
    int nDescandant = countVertexNumDown(node->C1) + countVertexNumDown(node->C2);
    *nVertex = nDescandant + nAncester;
    *separatorFrom = nDescandant;
    *vertexIdx = new unsigned int[*nVertex];
    unsigned int idxInUse = 0;
    SetUpVertexListDown(*vertexIdx, node->C1, idxInUse);
    SetUpVertexListDown(*vertexIdx, node->C2, idxInUse);
    SetUpVertexListUp(*vertexIdx, node->parent, idxInUse);  
}

void updateEdgeSuper(Supernode* node, std::mutex** mutexBlock, int sBlock, unsigned int** apsp){
    using namespace std;
    unsigned int* vertexIdx = nullptr;
    unsigned int* vertexIdxOwn = node->VertexIdx;
    int nVertex = 0;
    int separatorFrom = 0;
    int ownBlock = node->length;
    int n_b = (nVertex + sBlock - 1) / sBlock;
    SetUpVertexList(&vertexIdx, node, &nVertex, &separatorFrom);
    
    //stage 1
    for(int k = 0; k < ownBlock; k++){
        for(int i = 0; i < ownBlock; ++i){
            for(int j = 0; j < ownBlock; ++j){
                updateEdge(vertexIdxOwn[i],vertexIdxOwn[j],vertexIdxOwn[k],nVertex,apsp);
            }  
        }
    }
    //Stage 2
    {
        std::queue<std::thread> threads;
        for(int m = 0; m < n_b; ++m){
            unsigned int fromK = 0;
            unsigned int toK = ownBlock; 
            unsigned int fromI = m * sBlock;
            unsigned int toI = min((m + 1) * sBlock, nVertex); 
            unsigned int fromJ = fromK;
            unsigned int toJ = toK;
            threads.emplace(updateEdgeNoded, vertexIdx, vertexIdxOwn, vertexIdxOwn,
                fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp);
        }
        for(int m = 0; m < n_b; ++m){
            unsigned int fromK = 0;
            unsigned int toK = ownBlock;
            unsigned int fromI = fromK;
            unsigned int toI = toK;
            unsigned int fromJ = m * sBlock;
            unsigned int toJ = min((m + 1) * sBlock, nVertex); 
            threads.emplace(updateEdgeNoded, vertexIdxOwn, vertexIdx, vertexIdxOwn,
                fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp);
        }
        while(!threads.empty()){
            threads.front().join();
            threads.pop();
        }
    }
    //Stage 3
    {
        std::queue<std::thread> threads;
        for(int m = 0; m < n_b * n_b; ++m){
            int m_x = 0, m_y = 0;
            unsigned int range = getRange(m, n_b);
            unsigned int localIdx = m - getValue(range, n_b);
            m_x += range/4;
            m_y += range/4;
            unsigned int cL = n_b - (int)(range / 4) * 2;//current Length
            switch (range % 4){
            case 0:
                m_x += localIdx; break;
            case 1:
                m_x += (cL - 1); m_y += 1 + localIdx; break;
            case 2:
                m_x += (cL - 2) - localIdx; m_y += (cL - 1); break;
            case 3:
                m_y += (cL - 2) - localIdx; break;
            }
            unsigned int fromK = 0;
            unsigned int toK = ownBlock; 
            unsigned int fromI = m_x * sBlock;
            unsigned int toI = min((m_x + 1)*sBlock, nVertex); 
            unsigned int fromJ = m_y * sBlock;
            unsigned int toJ = min((m_y + 1)*sBlock, nVertex); 
            threads.emplace(updateEdgeNoded, vertexIdx, vertexIdx, vertexIdxOwn,
                fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp);
        }
        while(!threads.empty()){
            threads.front().join();
            threads.pop();
        }
    }
    
    delete[] vertexIdx;
}

Supernode* nestedDissectionBFS(unsigned int** distance, Supernode* node){
    //Return the root
    unsigned int* vertex = node->VertexIdx;
    unsigned int nVertexInUse = node->length;//number of vertecies in use
    if(nVertexInUse < 3){
        return node;
    }
    unsigned int nEdges[nVertexInUse];
    for(int i = 0; i < nVertexInUse; ++i){
        nEdges[i] = 0;
    }
    for(int i = 0; i < nVertexInUse; ++i){
        for(int j = 0; j < nVertexInUse; ++j){
            if((distance[vertex[i]][vertex[j]] != INF) && (i != j)){
                ++nEdges[i];
                ++nEdges[j];
            }
        }   
    }
    int minV = nEdges[0];
    for(int i = 1; i < nVertexInUse; ++i){
        if(minV > nEdges[i])
            minV = nEdges[i];
    }//smallest outgoing edge
    int candidateNum = 0;
    for(int i = 0; i < nVertexInUse; ++i){
        if(minV == nEdges[i])
            ++candidateNum;
    }//smallest outgoing edge
    
    int candidateMin[candidateNum];
    int candidateIdx = 0;
    for(int i = 0; i < nVertexInUse; ++i){
        if(minV == nEdges[i]){
            candidateMin[candidateIdx] = i;
            ++candidateIdx;
        }
    }//smallest outgoing edge
    
    //Ideality : (|C1|/|C2| + |C1|/|C2|) * |S|
    float ideality = MAXFLOAT;
    unsigned int idealDepthFrom = 0; //IdealDepth to use for seprator
    unsigned int idealDepthTo = 0; // Both inclusive
    unsigned int idealCandidate;
    unsigned int nS = 0, nC1 = 0, nC2 = 0;//actual optimal C1,C2,S
    unsigned int nDepth[nVertexInUse];
    unsigned int inUseVertex[nVertexInUse]; //Queue to check depth
    unsigned int separatorSize[nVertexInUse]; //If depth of idx chooses as a sepeartor, size of sepearator
    
    int minIdx;
    for(int candidateIdx = 0; candidateIdx < candidateNum; ++candidateIdx)
    {
        minIdx = candidateMin[candidateIdx];
        unsigned int depthed = 1; //How many elements have their depth in ndepth

        //This formula will be used to choose good dissection
        unsigned int usingIdx; //currently using Idx
        unsigned int writingIdx;//currently writing at
        for(int i = 0; i < nVertexInUse; ++i){
            separatorSize[i] = 0;
        }
        for(int i = 0; i < nVertexInUse; ++i){
            nDepth[i] = INF;
            inUseVertex[i] = -1;
        }
        nDepth[minIdx] = 0;
        inUseVertex[0] = minIdx;
        separatorSize[0] = 1;
        usingIdx = 0;
        writingIdx = 0;
        while(depthed < nVertexInUse){
            if(inUseVertex[usingIdx] == -1){
                break; //unreachable lefts
            }
            for(int j = 0; j < nVertexInUse; ++j){
                if(nDepth[j] == INF){
                    if(distance[vertex[inUseVertex[usingIdx]]][vertex[j]] != INF ||
                            distance[vertex[j]][vertex[inUseVertex[usingIdx]]] != INF){
                        nDepth[j] = nDepth[inUseVertex[usingIdx]] + 1;
                        ++separatorSize[nDepth[j]];
                        ++depthed;
                        inUseVertex[++writingIdx] = j;
                    }
                }
            }
            ++usingIdx;
        }
        unsigned int snC1 = separatorSize[0], snC2, snS; // just for computation
        for(int i = 1; i < nVertexInUse; ++i){
            if(separatorSize[i] == 0)
                break;
            snS = 0;
            for(int j = i; j < nVertexInUse; ++j){
                snS += separatorSize[j];    
                if(separatorSize[i] == 0)
                    break;
                snC2 = nVertexInUse - snC1 - snS;
                float currentIdeality = snS/float(snC1 + snC2) * abs(snC1/(float)snC2 + snC2/(float)snC1);
                if(ideality > currentIdeality){
                    idealCandidate = candidateMin[candidateIdx];
                    ideality = currentIdeality;
                    idealDepthFrom = i;
                    idealDepthTo = j;
                    nC1 = snC1;
                    nC2 = snC2;
                    nS = snS;
                }
            }
            snC1 += separatorSize[i];
        }
    }
    if(nC1 == 0 || nC2 == 0 || nS == 0){
        return node; //un dividable.
    }

    if(minIdx != idealCandidate){
        minIdx = idealCandidate;
        unsigned int depthed = 1; //How many elements have their depth in ndepth

        //This formula will be used to choose good dissection
        unsigned int usingIdx; //currently using Idx
        unsigned int writingIdx;//currently writing at
        unsigned int separatorSize[nVertexInUse]; //If depth of idx chooses as a sepeartor, size of sepearator
        for(int i = 0; i < nVertexInUse; ++i){
            separatorSize[i] = 0;
        }
        for(int i = 0; i < nVertexInUse; ++i){
            nDepth[i] = INF;
            inUseVertex[i] = -1;
        }
        nDepth[minIdx] = 0;
        inUseVertex[0] = minIdx;
        separatorSize[0] = 1;
        usingIdx = 0;
        writingIdx = 0;
        while(depthed < nVertexInUse){
            if(inUseVertex[usingIdx] == -1){
                break; //unreachable lefts
            }
            for(int j = 0; j < nVertexInUse; ++j){
                if(nDepth[j] == INF){
                    if(distance[vertex[inUseVertex[usingIdx]]][vertex[j]] != INF ||
                            distance[vertex[j]][vertex[inUseVertex[usingIdx]]] != INF){
                        nDepth[j] = nDepth[inUseVertex[usingIdx]] + 1;
                        ++separatorSize[nDepth[j]];
                        ++depthed;
                        inUseVertex[++writingIdx] = j;
                    }
                }
            }
            ++usingIdx;
        }
    }
    // Deviding a node to triple
    unsigned int newS[node->length];
    Supernode* C1 = new Supernode(SUPERNODE_ND_BFS);
    Supernode* C2 = new Supernode(SUPERNODE_ND_BFS);
    C1->length = nC1;
    C1->VertexIdx = node->VertexIdx + nS;
    C1->parent = node;
    C2->length = nC2;
    C2->VertexIdx = node->VertexIdx + nS + nC1;
    C2->parent = node;
    unsigned int idxC1 = nS,idxC2 = nS + nC1,idxS = 0;
    for(int i = 0; i < nVertexInUse; ++i){
        if(nDepth[i] < idealDepthFrom){
            newS[idxC1] = vertex[i];
            ++idxC1;
        }
        else if(nDepth[i] > idealDepthTo){
            newS[idxC2] = vertex[i];
            ++idxC2;
        }
        else{
            newS[idxS] = vertex[i];
            ++idxS;
        }
    }
    for(int i = 0; i < nVertexInUse; ++i){
        node->VertexIdx[i] = newS[i];
    }
    node->length = nS;
    node->C1 = C1;
    node->C2 = C2;
    std::thread t1 = std::thread(nestedDissectionBFS,distance,C1);
    nestedDissectionBFS(distance,C2);
    t1.join();
    return node;
}

Supernode* nestedDissectionCutMetis(unsigned int** distance, unsigned int** depthTable, Supernode* node, unsigned int from){
    // printf("callIn\n");
    unsigned int nVertex = node->length;
    if(nVertex < 3){ //Tosmall
        // printf("calledOut 1\n");
        return node;
    }
    unsigned int maxBox = 0;
    unsigned int maxI = 1, maxJ = 0;
    float ideality = MAXFLOAT; 
    for(int i = from + 1; i < from + nVertex - 1; ++i){
        for(int j = i; j < from + nVertex - 1; ++j){
            if(depthTable[i][j] == i){
                // unsigned int currentBox = (i - from) * (j - i + 1);
                unsigned int snC1 = i - from;
                unsigned int snC2 = j - i;
                unsigned int snS =  nVertex - snC1 - snC2;
                float currentIdeality = snS/float(snC1 + snC2) * abs(snC1/(float)snC2 + snC2/(float)snC1);
                if(ideality > currentIdeality){
                    ideality = currentIdeality;
                    maxI = i;
                    maxJ = j;
                }
            }
            else{
                break;
            }
        }  
    }
    if(maxI > maxJ){
        // printf("calledOut 2\n");
        return node;
    }//No possible saprsity exists
    //from(In) ~ maxI(Ex) C1
    //maxI(In) ~ maxJ(Ex) C2
    //maxJ(Ex) ~ from + nVertex (Ex) S

    unsigned int nC1 = maxI - from;
    unsigned int nC2 = maxJ - maxI;
    unsigned int nS =  nVertex - nC1 - nC2;
    Supernode* C1 = new Supernode(SUPERNODE_ND_MET);
    Supernode* C2 = new Supernode(SUPERNODE_ND_MET);
    C1->length = nC1;
    C1->VertexIdx = node->VertexIdx;
    C1->parent = node;

    C2->length = nC2;
    C2->VertexIdx = node->VertexIdx + nC1;
    C2->parent = node;
   
    node->VertexIdx = node->VertexIdx + nC1 + nC2;
    node->length = nS;
    node->C1 = C1;
    node->C2 = C2;

    // std::thread t1 = std::thread(nestedDissectionCutMetis, distance, depthTable, C1, from, threshold);
    nestedDissectionCutMetis(distance, depthTable, C1, from);
    nestedDissectionCutMetis(distance, depthTable, C2, from + nC1);
    
    // t1.join();

    // printf("calledOut 3\n");
    return node;
}

void DrawRedBox_real(unsigned char*** image, unsigned int from, unsigned int to,bool blue = false, unsigned int innerFrom = 0){
    if(blue){
        for(int i = from; i < innerFrom; ++i){
            image[i][innerFrom][2] = 255;
            image[i][innerFrom][1] = 0;
            image[i][innerFrom][0] = 0;
            image[innerFrom][i][2] = 255;
            image[innerFrom][i][1] = 0;
            image[innerFrom][i][0] = 0;
            image[i][to - 1][2] = 255;
            image[i][to - 1][1] = 0;
            image[i][to - 1][0] = 0;
            image[to - 1][i][2] = 255;
            image[to - 1][i][1] = 0;
            image[to - 1][i][0] = 0;
        }
        image[innerFrom][innerFrom][2] = 255;
        image[innerFrom][innerFrom][1] = 0;
        image[innerFrom][innerFrom][0] = 0;
        for(int i = innerFrom; i < to; ++i){
            image[from][i][2] = 255;
            image[from][i][1] = 0;
            image[from][i][0] = 0;
            image[i][from][2] = 255;
            image[i][from][1] = 0;
            image[i][from][0] = 0;
            image[i][to - 1][2] = 255;
            image[i][to - 1][1] = 0;
            image[i][to - 1][0] = 0;
            image[to - 1][i][2] = 255;
            image[to - 1][i][1] = 0;
            image[to - 1][i][0] = 0;
        }
    }
    else{
        for(int i = from; i < to; ++i){
            image[i][from][0] = 255;
            image[i][from][1] = 0;
            image[i][from][2] = 0;
            image[from][i][0] = 255;
            image[from][i][1] = 0;
            image[from][i][2] = 0;
            image[i][to - 1][0] = 255;
            image[i][to - 1][1] = 0;
            image[i][to - 1][2] = 0;
            image[to - 1][i][0] = 255;
            image[to - 1][i][1] = 0;
            image[to - 1][i][2] = 0;
        }
    }
}

void DrawRedBox(Supernode* node, unsigned char*** image, unsigned int to){
    if(node == nullptr )
        return;
    unsigned int nC1 = countEleNum(node->C1);
    unsigned int nC2 = countEleNum(node->C2);
    unsigned int nS = (node->length);
    if(node->C1 != nullptr && node->C2 != nullptr){
        DrawRedBox_real(image,to - nC1 - nC2 - nS, to, true, to - nS);
        DrawRedBox_real(image,to - nS, to);
        // DrawRedBox_real(image,to - nC2 - nS,to - nS);
        // DrawRedBox_real(image,to - nC1 - nC2 - nS,to - nC2 - nS);
        DrawRedBox(node->C1,image,to-nS-nC2);
        DrawRedBox(node->C2,image,to-nS);
    }
    else{
        DrawRedBox_real(image,to- nC1 - nC2 - nS, to);
    }
}

void produceVertexList(Supernode* node, std::vector<unsigned int>& list){
    if(node == nullptr)
        return;
    produceVertexList(node->C1, list);
    produceVertexList(node->C2, list);
    for(int i = 0; i < node->length; ++i)
        list.emplace_back(node->VertexIdx[i]);
}

void OutPutMatrixFile(const char* fileName, unsigned int nVertex, unsigned int** distance, Supernode* root){
    // === Matrix printing ===//
    unsigned char*** resultIamge;
    unsigned int* vertex;
    std::vector<unsigned int> vertexList;
    if(root != nullptr){
        produceVertexList(root,vertexList);
        if(vertexList.size() != nVertex){
            printf("wrong output Matrix Image\n");
            exit(0);
        }
    }
    vertex = new unsigned int[nVertex];
    if(root == nullptr){
        for(int i = 0; i < nVertex; ++i){
            vertex[i] = i;
        }
    }
    else{
        for(int i = 0; i < nVertex; ++i){
            vertex[i] = vertexList[i];
        }
    }
    resultIamge = new unsigned char**[nVertex];
    for(int i = 0; i<nVertex; ++i){
        resultIamge[i] = new unsigned char*[nVertex];
        for(int j = 0; j<nVertex; ++j){
            resultIamge[i][j] = new unsigned char[3];
        }
    }
    for(int i = 0; i<nVertex; ++i){
        for(int j = 0; j<nVertex; ++j){
            if(distance[vertex[i]][vertex[j]] == INF){
                resultIamge[i][j][0] = 255;
                resultIamge[i][j][1] = 255;
                resultIamge[i][j][2] = 255;
            }
            else{
                resultIamge[i][j][0] = 0;
                resultIamge[i][j][1] = 0;
                resultIamge[i][j][2] = 0;
            }
        }
    }
    if(root != nullptr)
        DrawRedBox(root,resultIamge,nVertex);
    std::ofstream matrixImage(fileName);
    matrixImage.write("P3\n",3);
    std::string headerMat = "";
    headerMat += std::to_string(nVertex) + " " 
                + std::to_string(nVertex) + "\n" ;
    matrixImage.write(headerMat.data(), headerMat.size());
    matrixImage.write("255\n",4);
    for(int i = 0; i<nVertex; ++i){
        for(int j = 0; j<nVertex; ++j){
            if(j != 0){
                matrixImage.write("  ",2);
            }
            std::string tempString = std::to_string(resultIamge[i][j][0]) + " "
                                    + std::to_string(resultIamge[i][j][1]) + " "
                                    + std::to_string(resultIamge[i][j][2]);
            matrixImage.write(tempString.data(), tempString.size());
        }matrixImage.write("\n", 1);
    }
    matrixImage.close();
    for(int i = 0; i < nVertex; ++i){
        for(int j = 0; j < nVertex; ++j){
            delete[] resultIamge[i][j];
        }
        delete[] resultIamge[i];
    }
    delete[] resultIamge;
    delete[] vertex;
    // === Matrix printing ===//
}

Supernode* nestedDissectionMetis(unsigned int** distance, Supernode* root){
    
    unsigned long long from = PAPI_get_real_usec();
    idx_t nVertex = (idx_t)(root->length);
    idx_t* rowPtr = new idx_t[nVertex + 1];
    std::vector<idx_t> nonZero;
    idx_t nEdges = 0;
    for(idx_t i = 0; i < nVertex; ++i){
        rowPtr[i] = nEdges;
        for(idx_t j = 0; j < nVertex; ++j){
            if(((distance[i][j] != INF) || (distance[j][i] != INF)) && (i != j)){
                nonZero.emplace_back(j);
                ++nEdges;
            }
        }
    }
    rowPtr[nVertex] = nEdges;
    idx_t* colIdx = new idx_t[nEdges];
    for(int i = 0; i < nEdges; ++i){
        colIdx[i] = nonZero[i];
    }
    idx_t* perm  = new idx_t[nVertex];
    idx_t* iperm = new idx_t[nVertex];
    
    
    unsigned long long to = PAPI_get_real_usec();
    printf("Post took : %lld\n", to - from);

    // printf("n : %lld m : %lld\n", nVertex, nEdges);
    // for(int i = 0; i <= nVertex; ++i){
    //     printf("%llu ", rowPtr[i]);
    // }printf("\n");

    // for(int i = 0; i < nEdges; ++i){
    //     printf("%llu ", colIdx[i]);
    // }printf("\n");

    from = PAPI_get_real_usec();
    int resultND;
    resultND = METIS_NodeND(&nVertex,rowPtr,colIdx,NULL,NULL,perm,iperm);

    if(resultND != METIS_OK){
        printf("ERROR CODE : %d\n",resultND);
        exit(0);
        return root;
    }
    to = PAPI_get_real_usec();
    printf("ND took : %lld\n", to - from);
    
    from = PAPI_get_real_usec();
    unsigned int* vertex = root->VertexIdx;
    for(int i = 0; i < nVertex; ++i){
        vertex[i] = perm[i];
    }
    
    unsigned int depth[nVertex];
    for(int i = 0; i < nVertex; ++i){
        depth[i] = 0;
        for(int j = 0; j < i; ++j){
            if(distance[vertex[i]][vertex[j]] != INF || distance[vertex[j]][vertex[i]] != INF){
                break;
            }
            ++depth[i];
        }
    }


    unsigned int* depthTable[nVertex] = {0};
    for(int i = 0; i < nVertex; ++i){
        depthTable[i] = new unsigned int[nVertex];
        for(int j = 0; j < nVertex; ++j){
            depthTable[i][j] = nVertex;
        }
    }
    
    for(int i = 0; i < nVertex; ++i){
        depthTable[i][i] = depth[i];
    }

    for(int i = 1; i < nVertex; ++i){
        for(int j = 0; j < nVertex - i; ++j){
            depthTable[j][j + i] = std::min(depthTable[j][j + i - 1], depthTable[j + 1][j + i]);
        }
    }//The smallest depth

    // printf("===Distance===\n");
    // for(int i = 0; i < nVertex; ++i){
    //     for(int j = 0; j < nVertex; ++j){
    //         printf("%c ",distance[vertex[i]][vertex[j]] == INF? '0' : '1');
    //     }printf("\n");
    // }

    // printf("===Depth===\n");
    // for(int i = 0; i < nVertex - 1; ++i){
    //     printf("%llu ",depth[i]);
    // }printf("\n");
    
    // printf("===DepthTable===\n");
    // for(int i = 0; i < nVertex - 1; ++i){
    //     for(int j = 0; j < nVertex - 1; ++j){
    //         printf("%llu ",depthTable[i][j]);
    //     }printf("\n");
    // }
    
    // printf("====BoxTable====\n");
    // for(int i = 0; i < nVertex - 1; ++i){
    //     for(int j = 0; j < nVertex - 1; ++j){
    //         printf("%llu ",depthTable[i][j] * (j - i + 1));
    //     }printf("\n");
    // }

    //depthTable[i][j] * (j - i + 1) means how big saprsity exists at matrix.
    // i <= j needed.

    nestedDissectionCutMetis(distance, depthTable, root, 0);
    
    

    for(int i = 0; i < nVertex; ++i){
        delete[] depthTable[i];
    }
    delete[] rowPtr;
    delete[] colIdx;
    delete[] perm;
    delete[] iperm;
    to = PAPI_get_real_usec();
    printf("ND cut took : %lld\n", to - from);
    from = PAPI_get_real_usec();
    
    return root;
}

Supernode* nestedDissection(unsigned int** distance, int nVertex, unsigned char MODE){
    Supernode* result = new Supernode(MODE);
    result->C1 = nullptr;
    result->C2 = nullptr;
    result->VertexIdx = new unsigned int[nVertex];
    for(int i = 0; i < nVertex; ++i)
      result->VertexIdx[i] = i;
    result->length = nVertex;
    switch (MODE)
    {
    case SUPERNODE_ND_MET:
        nestedDissectionMetis(distance, result);
        break;
    case SUPERNODE_ND_BFS:
        nestedDissectionBFS(distance, result);
        break;
    }
    return result;
}

void PushToScheduler(Supernode* root, std::stack<std::pair<Supernode*,unsigned int>>& operationStack){
    std::queue<std::pair<Supernode*,unsigned  int>> bfsqueue;
    bfsqueue.push(std::pair<Supernode*, unsigned int>(root, 0));
    while(!bfsqueue.empty()){
        std::pair<Supernode*,unsigned int> top = bfsqueue.front();
        operationStack.emplace(top);
        if(top.first->C1 != nullptr)
            bfsqueue.push(std::pair<Supernode *, unsigned int>(top.first->C1, top.second + 1));
        if(top.first->C2 != nullptr)
            bfsqueue.push(std::pair<Supernode *, unsigned int>(top.first->C2, top.second + 1));
        bfsqueue.pop();
    }
}

char readByChar(){
    std::string input;
    std::cin >> input;
    return input[0];
}

unsigned int readByUI(){
    std::string input;
    std::cin >> input;
    unsigned int ret = 0;
    for(int i = 0; i < input.length(); ++i){
        char ch = input[i];
        if('0' <= ch && ch <= '9'){
            ret *= 10;
            ret += ch - '0';
        }
        else{
            return 0;
        }
    }
    return ret;
}

unsigned int readByBinary(){
    std::string input;
    std::cin >> input;
    unsigned int ret = 0;
    for(int i = 0; i < input.length(); ++i){
        char ch = input[i];
        if('0' <= ch && ch <= '1'){
            ret *= 2;
            ret += ch - '0';
        }
        else{
            return 0;
        }
    }
    return ret;
}

int main(){
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if(retval != PAPI_VER_CURRENT) exit(0);
    srand(time(NULL));
    using namespace std;
    unsigned int nVertex = 0, nEdges = 0,sBlock = 0,nThreads = 0,superFWthreshold_BFS = 0, superFWthreshold_MET = 0, option;
    unsigned int graphType;
    unsigned int** distance = nullptr;
    char terminate = 0;
    bool skipClassic;
    //Criteria = blocked threading
    // globalThreadCall[0].store(0);
    // globalThreadCall[1].store(0);
    // globalThreadCall[2].store(0);
    vector<int> speedUpResultsBlocked;
    vector<int> speedUpResultsPerLine;
    vector<int> speedUpResultsPerElem;
    vector<int> speedUpResultsPerThread;
    vector<int> speedUpResultsBlockedqMul;
    vector<int> speedUpResultsSuperFWBFS;
    vector<int> speedUpResultsSuperFWMET;
    vector<int> speedUpResultsPathDoubling;
    int trial = 0;
    unsigned long long classicTime;
    unsigned long long timeConsume;
    constexpr unsigned int BlockedFlag = 128;
    constexpr unsigned int PerLineFlag = 64;
    constexpr unsigned int PerElemFlag = 32;
    constexpr unsigned int PerThreFlag = 16;
    constexpr unsigned int MulBlocFlag = 8;
    constexpr unsigned int SuperFWFlagBFS = 4;
    constexpr unsigned int SuperFWFlagMET = 2;
    constexpr unsigned int PathDubFlag = 1;
    std::string graphName;
    do{
        // Make a new graph input
        if(terminate != 'r'){
            unsigned int approxiGap = 0;
            speedUpResultsBlocked.clear();
            speedUpResultsPerLine.clear();
            speedUpResultsPerElem.clear();
            speedUpResultsPerThread.clear();
            speedUpResultsBlockedqMul.clear();
            speedUpResultsSuperFWBFS.clear();
            speedUpResultsSuperFWMET.clear();
            speedUpResultsPathDoubling.clear();
            ifstream graphFile;
            string fileName;
            if(distance != nullptr){   
                for(int i = 0; i < nVertex; ++i){
                    delete[] distance[i];
                }
                delete[] distance;
            }
            while(true){
                approxiGap = 0;
                cout << "Testset classic(fixed) blocked perLine perElement perThread blockThread SuperFWBFS superFWMET Pathdoubling)"<< endl;
                cout << "eg1) blocked perline perElement : 11100000, eg2)superFMETWonly : 10"<< endl;
                cout <<"Binary format (0 to exit) : ";
                option = readByBinary();
                if(option == 0)
                    exit(0);
                cout << "Graph type " << endl;
                cout << "0) linear, 1)EB, 2)app-linear, 3)tree, 4)dense, 5)true-random) 6)planer 7)Input file : ";
                graphType = readByUI();
                if(graphType == 7){
                    cout << "Type file name : ";
                    cin >> fileName;
                }
                if(graphType == 2){
                    cout << "Approximiation gap: ";
                    approxiGap = readByUI();
                }
                if(graphType == 6){
                    cout << "Size of the plane: ";
                    nVertex = readByUI();
                    nVertex = nVertex*nVertex;
                }
                else if(!(graphType == 7)){
                    cout << "Number of vertices: ";
                    nVertex = readByUI();
                }
                if(!(graphType == 0 || graphType == 3 || graphType == 4 || graphType == 6 || graphType == 7)){
                    cout << "Number of edges: ";
                    nEdges = readByUI();
                }
                if((BlockedFlag|MulBlocFlag)&option){
                    cout << "Size of block: ";
                    sBlock = readByUI();
                }
                if(PerThreFlag&option){
                    cout << "Number of threads: ";
                    nThreads = readByUI();
                }
                if(SuperFWFlagBFS&option){
                    cout << "Size of threshold of superFW(BFS): ";
                    superFWthreshold_BFS = readByUI();
                }
                if(SuperFWFlagMET&option){
                    cout << "Size of threshold of superFW(MET): ";
                    superFWthreshold_MET = readByUI();
                }
                char tempCharIn;
                
                if(option&MulBlocFlag){
                    do{
                        cout << "Skip untill Multithreaded Block Fw (y or n): ";
                        tempCharIn = readByChar();
                    }while(tempCharIn != 'y' && tempCharIn != 'n');
                    skipClassic = tempCharIn == 'y';
                }
                else{
                    skipClassic = false;
                }
                if(skipClassic){
                    option |= MulBlocFlag;
                    if(BlockedFlag & option)
                       option -= BlockedFlag;
                    if(PerElemFlag & option)
                       option -= PerElemFlag;
                    if(PerLineFlag & option)
                       option -= PerLineFlag;
                    if(PerThreFlag & option)
                       option -= PerThreFlag;
                }
                
                bool tryAgain = false;
                if(graphType == 7){
                    cout << ("Try to load ../Testcases/" + fileName + ".mtx") << endl;
                    graphFile.open("../Testcases/" + fileName + ".mtx");
                    if(graphFile.is_open()){
                        string line;
                        while(getline(graphFile,line)){
                            if(line[0] != '%')
                                break;
                        }
                        nVertex = 0;
                        int nVertex2 = 0;
                        int indexOfLine = 0;
                        for(; indexOfLine < line.length(); ++indexOfLine){
                            if(line[indexOfLine] == ' '){
                                ++indexOfLine;
                                break;
                            }
                            nVertex *= 10;
                            nVertex += line[indexOfLine] - '0';
                        }
                        for(;indexOfLine < line.length(); ++indexOfLine){
                            if(line[indexOfLine] == ' ' || line[indexOfLine] == '\n'){
                                ++indexOfLine;
                                break;
                            }
                            nVertex2 *= 10;
                            nVertex2 += line[indexOfLine] - '0';
                        }
                        if(nVertex2 != nVertex)
                            tryAgain = true;
                    }
                    else{
                        cout << "failed to open" << endl;
                        tryAgain = true;
                    }
                }
                
                if(tryAgain){
                    cout << "Invalid file input" << endl;
                }
                else if(nEdges == 0 && (!(graphType == 0 || graphType == 3 || graphType == 4 || graphType == 6 || graphType == 7 )))
                    cout << "Please make valid graph" << endl;
                else if(nVertex == 0 && graphType != 7)
                    cout << "Please make valid graph" << endl;
                else if((nEdges < nVertex) && (!(graphType == 0 || graphType == 3 || graphType == 4 || graphType == 5 || graphType == 6 || graphType == 7)))
                    cout << "Please type more edges than vertex" << endl;
                else if(nEdges > nVertex*nVertex && (!(graphType == 0 || graphType == 3 || graphType == 4|| graphType == 6 || graphType == 7)))
                    cout << "Please type valid number of edges" << endl;
                else if((sBlock > nVertex) && ((BlockedFlag|MulBlocFlag) & option))
                    cout << "The block size(" << sBlock << ") should be smaller than or queal to the number of vertex (" << nVertex << ")" << endl;
                else if((nThreads < 1) && (PerThreFlag & option))
                    cout << "At least one thread required" << endl;
                else
                    break;
            }
            INF = nVertex * 2;
            distance = new unsigned int*[nVertex];
            for(int i = 0; i < nVertex; ++i){
                distance[i] = new unsigned int[nVertex];
                for(int j = 0; j < nVertex; ++j){
                    distance[i][j] = INF;
                }
            }
            printf("Generating graph\n");
            switch(graphType){
            case 0: // Linear Graph
                nEdges = 0;
                graphName = "Linear graph";
                for(int i = 0; i < nVertex - 1; ++i){
                    nEdges += (distance[i][i+1] == INF);
                    nEdges += (distance[i+1][i] == INF);
                    distance[i][i+1] = 1;
                    distance[i+1][i] = 1;
                }
                break;
            case 1: // BA Graph
                graphName = "Barbasi-Albert graph";
                {
                    unsigned int m = nEdges;
                    for(int i = 1; i < nVertex; ++i,--m){
                        // connectivity guarantee
                        // i will be connected with to
                        while(true){
                            unsigned int to = rand() % i;
                            if(to != i && distance[i][to] == INF){
                                distance[i][to] = 1;
                                distance[to][i] = 1;
                                break;
                            }
                        }
                    }
                    for(;m>0;--m){
                        while(true){
                            unsigned int from = rand() % nVertex;
                            unsigned int to = rand() % nVertex;
                            if(to != from && distance[from][to] == INF){
                                distance[from][to] = 1;
                                distance[to][from] = 1;
                                break;
                            }
                        } 
                    }
                }
                break;
            case 2: // Approximately linear Graph
                {
                    graphName = "Appriximately linear graph";
                    unsigned int m = nEdges - (nVertex - 1);
                    for(int i = 0; i < nVertex - 1; ++i){
                        distance[i][i+1] = 1;
                        distance[i+1][i] = 1;
                    }
                    for(;m>0;--m){
                        while(true){
                            unsigned int from = rand() % nVertex;
                            unsigned int to = from + (rand() % approxiGap - approxiGap / 2);
                            if(to < 0)
                                to = 0;
                            if(to >= nVertex)
                                to = nVertex - 1;
                            if(to != from && distance[from][to] == INF){
                                distance[from][to] = 1;
                                distance[to][from] = 1;
                                break;
                            }
                        } 
                    }
                }
                break;
            case 3: // Tree Graph
                graphName = "Tree graph";
                nEdges = 0;
                for(int i = 0; i < nVertex ; ++i){
                    if( 2 * i + 1 < nVertex){
                        nEdges += (distance[i][2 * i + 1] == INF);
                        distance[i][2 * i + 1] = 1;
                        nEdges += (distance[2 * i + 1][i] == INF);
                        distance[2 * i + 1][i] = 1;
                    }
                    if( 2 * i + 2 < nVertex){
                        nEdges += (distance[i][2 * i + 2] == INF);
                        distance[i][2 * i + 2] = 1;
                        nEdges += (distance[2 * i + 2][i] == INF);
                        distance[2 * i + 2][i] = 1;
                    }
                }
                break;
            case 4: // Dense Graph
                graphName = "Dense graph";
                nEdges = 0;
                for(int i = 0; i < nVertex ; ++i){
                    for(int j = 0; j < nVertex; ++j){
                        nEdges += (distance[i][j] == INF);
                        nEdges += (distance[j][i] == INF);
                        distance[i][j] = 1;
                        distance[j][i] = 1;
                    }
                }
                break;
            case 5: // True random Graph
                graphName = "True random graph";
                {
                    unsigned int m = nEdges;
                    for(;m>0;--m){
                        while(true){
                            unsigned int from = rand() % nVertex;
                            unsigned int to = rand() % nVertex;
                            if(to != from && distance[from][to] == INF){
                                distance[from][to] = 1;
                                distance[to][from] = 1;
                                break;
                            }
                        } 
                    }
                }
                break;
            case 6: // Linear Graph
                {
                    nEdges = 0;
                    graphName = "Planer graph";
                    unsigned int size = std::sqrt(nVertex);
                    for(int i = 0; i < nVertex; ++i){
                        if(i + size < nVertex){
                            nEdges += (distance[i][i+size] == INF);
                            nEdges += (distance[i+size][i] == INF);
                            distance[i][i+size] = 1;
                            distance[i+size][i] = 1;
                        }
                        if(i + 1 < nVertex){
                            nEdges += (distance[i][i+1] == INF);
                            nEdges += (distance[i+1][i] == INF);
                            distance[i][i+1] = 1;
                            distance[i+1][i] = 1;
                        }
                    }
                }
                break;
            case 7: // Given input
                {
                    nEdges = 0;
                    graphName = "Given graph : " + fileName;
                    string line;
                    while(getline(graphFile, line)){
                        int indexOfLine = 0;
                        int a, b;
                        a = 0;
                        b = 0;
                        for(; indexOfLine < line.length(); ++indexOfLine){
                            if(line[indexOfLine] == ' ' || line[indexOfLine] == '\n' ){
                                ++indexOfLine;
                                break;
                            }
                            a *= 10;
                            a += line[indexOfLine] - '0';
                        }
                        for(; indexOfLine < line.length(); ++indexOfLine){
                            if(line[indexOfLine] == ' ' || line[indexOfLine] == '\n' ){
                                ++indexOfLine;
                                break;
                            }
                            b *= 10;
                            b += line[indexOfLine] - '0';
                        }
                        a -= 1;
                        b -= 1;
                        nEdges += (distance[a][b] == INF);
                        distance[a][b] = 1;
                        nEdges += (distance[b][a] == INF);
                        distance[b][a] = 1;
                    }
                }
                graphFile.close();
                break;
            }
            // OutPutMatrixFile("./MatrixImage_Input.ppm",nVertex,distance,nullptr);
            printf("Graph generated\n");
        }
        unsigned int** apsp1;
        if(!skipClassic){
            initResult(nVertex,distance,&apsp1);
            printf("Memory prepared for classic\n");
        }

        unsigned int** apsp2;
        if(BlockedFlag&option){
            initResult(nVertex,distance,&apsp2);
            printf("Memory prepared for blocked sequential\n");
        }

        unsigned int** apsp3;
        if(PerLineFlag&option){
            initResult(nVertex,distance,&apsp3);
            printf("Memory prepared for threading per line\n");
        }

        unsigned int** apsp4;
        if(PerElemFlag&option){
            initResult(nVertex,distance,&apsp4);
            printf("Memory prepared for threading per element\n");
        }

        unsigned int** apsp5;    
        if(PerThreFlag&option){
            initResult(nVertex,distance,&apsp5);
            printf("Memory prepared for threading per thread\n");
        }

        unsigned int** apsp6;  
        if(MulBlocFlag&option){
            initResult(nVertex,distance,&apsp6);
            printf("Memory prepared for Blocked threading\n");
        }

        unsigned int** apsp7;  
        if(SuperFWFlagBFS&option){
            initResult(nVertex,distance,&apsp7);
            printf("Memory prepared for SuperFWBFS\n");
        }

        unsigned int** apsp9;          
        if(SuperFWFlagMET&option){
            initResult(nVertex,distance,&apsp9);
            printf("Memory prepared for SuperFWMET\n");
        }
    
        unsigned int* apsp8;  
        if(PathDubFlag&option){
            initResultGPU(nVertex,distance,&apsp8);
            printf("Memory prepared for Path doubling\n");
        }

        //classical sequential
        if(!skipClassic)
        {
            unsigned long long from = PAPI_get_real_usec();
            for(int k = 0; k < nVertex; ++k){
                for(int i = 0; i < nVertex; ++i){
                    for(int j = 0; j < nVertex; ++j){
                        updateEdge(i,j,k,nVertex,apsp1);
                    }  
                }  
            }
            unsigned long long to = PAPI_get_real_usec();
            classicTime = (to-from);
            cout << "classic : " << classicTime << "usec" << endl;
        }

        //blocked sequential
        if((!skipClassic) && (BlockedFlag&option))
        {
            unsigned long long from = PAPI_get_real_usec();
            unsigned int n_b = (nVertex + sBlock - 1)/sBlock;
            for(int b = 0; b < n_b; ++b){
                //stage 1
                unsigned int fromK = b * sBlock;
                unsigned int toK = min((b + 1) * sBlock, nVertex);
                unsigned int fromI = fromK;
                unsigned int toI = toK;
                unsigned int fromJ = fromK;
                unsigned int toJ = toK;
                for(int k = fromK; k < toK; k++){
                    for(int i = fromI; i < toI; ++i){
                        for(int j = fromJ; j < toJ; ++j){
                            updateEdge(i,j,k,nVertex,apsp2);
                        }  
                    }
                }
                //Stage 2
                for(int m = 0; m < n_b; ++m){
                    if(m == b)
                        continue;
                    unsigned int fromK = b * sBlock;
                    unsigned int toK = min((b + 1) * sBlock, nVertex); 
                    unsigned int fromI = m * sBlock;
                    unsigned int toI = min((m + 1) * sBlock, nVertex); 
                    unsigned int fromJ = fromK;
                    unsigned int toJ = toK;
                    for(int k = fromK; k < toK; k++){
                        for(int i = fromI; i < toI; ++i){
                            for(int j = fromJ; j < toJ; ++j){
                                updateEdge(i,j,k,nVertex,apsp2);
                            }  
                        }
                    }
                }
                for(int m = 0; m < n_b; ++m){
                    if(m == b)
                        continue;
                    unsigned int fromK = b * sBlock;
                    unsigned int toK = min((b + 1) * sBlock, nVertex);
                    unsigned int fromI = fromK;
                    unsigned int toI = toK;
                    unsigned int fromJ = m * sBlock;
                    unsigned int toJ = min((m + 1) * sBlock, nVertex); 
                    for(int k = fromK; k < toK; k++){
                        for(int i = fromI; i < toI; ++i){
                            for(int j = fromJ; j < toJ; ++j){
                                updateEdge(i,j,k,nVertex,apsp2);
                            }  
                        }
                    }
                }
                //Stage 3
                for(int m = 0; m < n_b * n_b; ++m){
                    int m_x = 0, m_y = 0;
                    unsigned int range = getRange(m, n_b);
                    unsigned int localIdx = m - getValue(range, n_b);
                    m_x += range/4;
                    m_y += range/4;
                    unsigned int cL = n_b - (int)(range / 4) * 2;//current Length
                    switch (range % 4){
                    case 0:
                        m_x += localIdx; break;
                    case 1:
                        m_x += (cL - 1); m_y += 1 + localIdx; break;
                    case 2:
                        m_x += (cL - 2) - localIdx; m_y += (cL - 1); break;
                    case 3:
                        m_y += (cL - 2) - localIdx; break;
                    }
                    if(m_x == b)
                        continue;
                    if(m_y == b)
                        continue;
                    unsigned int fromK = b*sBlock;
                    unsigned int toK = min((b + 1)*sBlock, nVertex); 
                    unsigned int fromI = m_x * sBlock;
                    unsigned int toI = min((m_x + 1)*sBlock, nVertex); 
                    unsigned int fromJ = m_y * sBlock;
                    unsigned int toJ = min((m_y + 1)*sBlock, nVertex); 
                    for(int k = fromK; k < toK; k++){
                        for(int i = fromI; i < toI; ++i){
                            for(int j = fromJ; j < toJ; ++j){
                                updateEdge(i,j,k,nVertex,apsp2);
                            }  
                        }
                    }
                }
            }
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "Blocked : " <<timeConsume << "usec" << endl;
            speedUpResultsBlocked.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsBlocked.begin(),speedUpResultsBlocked.end());
        }
        
        //classical thread per line   
        if(PerLineFlag&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            for(int k = 0; k < nVertex; ++k){
                thread threads[nVertex];
                for(int j = 0; j < nVertex; ++j){
                    threads[j] = thread(updateEdgePerLine,j,k,nVertex,apsp3);
                }
                for(int j = 0; j < nVertex; ++j){
                    threads[j].join();
                }
            }
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "threading per line : " << timeConsume << "usec" << endl;
            speedUpResultsPerLine.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsPerLine.begin(),speedUpResultsPerLine.end());
        }

        //classical thread per  
        if(PerElemFlag&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            for(int k = 0; k < nVertex; ++k){
                thread threads[nVertex*nVertex];
                for(int i = 0; i < nVertex; ++i){
                    for(int j = 0; j < nVertex; ++j){
                        threads[i + j * nVertex] = thread(updateEdgePerElem,i,j,k,nVertex,apsp4);
                    }
                }
                for(int i = 0; i < nVertex; ++i){
                    for(int j = 0; j < nVertex; ++j){
                        threads[i + j * nVertex].join();
                    }
                }
            }
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "threading per elem : " << timeConsume << "usec" << endl;
            speedUpResultsPerElem.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsPerElem.begin(),speedUpResultsPerElem.end());
        }
    
        //classical thread per threads
        if(PerThreFlag&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            for(int k = 0; k < nVertex; ++k){
                thread threads[nVertex*nVertex];
                for(int i = 0; i < nThreads; ++i){
                    threads[i] = thread(updateEdgePerThread,i,k,nThreads,nVertex,apsp5);
                }
                for(int i = 0; i < nThreads; ++i){
                    threads[i].join();
                }
            }
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "threading per thread : " << timeConsume << "usec" << endl;
            speedUpResultsPerThread.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsPerThread.begin(),speedUpResultsPerThread.end());
        }
        
        //blocked threading
        if(MulBlocFlag&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            unsigned int n_b = (nVertex + sBlock - 1)/sBlock;
            for(int b = 0; b < n_b; ++b){
                //stage 1
                unsigned int fromK = b * sBlock;
                unsigned int toK = min((b + 1) * sBlock, nVertex);
                unsigned int fromI = fromK;
                unsigned int toI = toK;
                unsigned int fromJ = fromK;
                unsigned int toJ = toK;
                for(int k = fromK; k < toK; k++){
                    for(int i = fromI; i < toI; ++i){
                        for(int j = fromJ; j < toJ; ++j){
                            updateEdge(i,j,k,nVertex,apsp6);
                        }  
                    }
                }
                //Stage 2
                {
                    int threadCount = 0;
                    thread threads[n_b * 2 - 2];
                    for(int m = 0; m < n_b; ++m){
                        if(m == b)
                            continue;
                        unsigned int fromK = b * sBlock;
                        unsigned int toK = min((b + 1) * sBlock, nVertex); 
                        unsigned int fromI = m * sBlock;
                        unsigned int toI = min((m + 1) * sBlock, nVertex); 
                        unsigned int fromJ = fromK;
                        unsigned int toJ = toK;
                        threads[threadCount++] = thread(updateEdgeBlocked,fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp6);
                    }
                    for(int m = 0; m < n_b; ++m){
                        if(m == b)
                            continue;
                        unsigned int fromK = b * sBlock;
                        unsigned int toK = min((b + 1) * sBlock, nVertex);
                        unsigned int fromI = fromK;
                        unsigned int toI = toK;
                        unsigned int fromJ = m * sBlock;
                        unsigned int toJ = min((m + 1) * sBlock, nVertex); 
                        threads[threadCount++] = thread(updateEdgeBlocked,fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp6);
                    }
                    for(int i = 0; i < 2* n_b - 2; ++i){
                        threads[i].join();
                    }
                }
                //Stage 3
                {
                    int threadCount = 0;
                    thread threads[(n_b -1) * (n_b -1)];
                    for(int m = 0; m < n_b * n_b; ++m){
                        int m_x = 0, m_y = 0;
                        unsigned int range = getRange(m, n_b);
                        unsigned int localIdx = m - getValue(range, n_b);
                        m_x += range/4;
                        m_y += range/4;
                        unsigned int cL = n_b - (int)(range / 4) * 2;//current Length
                        switch (range % 4){
                        case 0:
                            m_x += localIdx; break;
                        case 1:
                            m_x += (cL - 1); m_y += 1 + localIdx; break;
                        case 2:
                            m_x += (cL - 2) - localIdx; m_y += (cL - 1); break;
                        case 3:
                            m_y += (cL - 2) - localIdx; break;
                        }
                        if(m_x == b)
                            continue;
                        if(m_y == b)
                            continue;
                        unsigned int fromK = b*sBlock;
                        unsigned int toK = min((b + 1)*sBlock, nVertex); 
                        unsigned int fromI = m_x * sBlock;
                        unsigned int toI = min((m_x + 1)*sBlock, nVertex); 
                        unsigned int fromJ = m_y * sBlock;
                        unsigned int toJ = min((m_y + 1)*sBlock, nVertex); 
                        threads[threadCount++] = thread(updateEdgeBlocked,fromI,toI,fromJ,toJ,fromK,toK,nVertex,apsp6);
                    }
                    for(int i = 0; i < (n_b - 1) * (n_b - 1); ++i){
                        threads[i].join();
                    }
                }
            }
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "Blocked threading : " <<timeConsume << "usec" << endl;
             if(skipClassic){
                classicTime = timeConsume;
            }
            else{
                speedUpResultsBlockedqMul.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
                sort(speedUpResultsBlockedqMul.begin(),speedUpResultsBlockedqMul.end());
            }
        }
        
        //SuperFWBFS
        unsigned int superFWdepth_BFS = 0;
        unsigned int superFWnum_BFS = 0;
        if(SuperFWFlagBFS&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            std::mutex** mutexBlock = nullptr;
            mutexBlock = new std::mutex*[nVertex];
            for(int i = 0; i < nVertex; ++i){
                mutexBlock[i] = new std::mutex[nVertex];
            }
            Supernode* root = nestedDissection(distance,nVertex,SUPERNODE_ND_BFS);
                
            stack<pair<Supernode*,unsigned int>> operationSchedule;
            PushToScheduler(root,operationSchedule); 
            superFWdepth_BFS = countDepth(root);
            superFWnum_BFS = countNum(root);
            //unsigned long long timeConsumePerDepth[superFWdepth] = {0,};
            //unsigned long long threadCalled[superFWdepth][3];
            //for(int i = 0; i < superFWdepth; ++i){
            //    threadCalled[i][0] = 0;
            //    threadCalled[i][1] = 0;
            //    threadCalled[i][2] = 0;
            //}
            //printBFS(root);
            unsigned long long toPre = PAPI_get_real_usec();
            while(!operationSchedule.empty()){
                //Schedule  =>  Buffer1   (Stage1)
                //Buffer1   =>  Buffer2 
                //Buffer2   =>  free      (Stage3)
                unsigned int currentDepth = operationSchedule.top().second;
                queue<std::thread> threadStack;
                while(true){//SuperRelax
                    threadStack.emplace(updateEdgeSuper, operationSchedule.top().first,mutexBlock, superFWthreshold_BFS,apsp7);
                    operationSchedule.pop();
                    if(operationSchedule.empty())
                        break;
                    if(operationSchedule.top().second != currentDepth)
                        break;
                }
                while(!threadStack.empty()){
                    threadStack.front().join();
                    threadStack.pop();
                }
            }
            //printBFS(root);
            // OutPutMatrixFile("./MatrixImage_BFS.ppm",nVertex,distance,root);
            delete root;
            for(int i = 0; i < nVertex; ++i){
                    delete[] mutexBlock[i];
                }
            delete[] mutexBlock;
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "SuperFW threading BFS (Post :  " << double(toPre - from)/(to - from) * 100 <<"%) : " << timeConsume << "usec" << endl;
            
            // for(int i = 0; i < superFWdepth-1; ++i){
            //    threadCalled[i][0] -= threadCalled[i+1][0];
            //    threadCalled[i][1] -= threadCalled[i+1][1];
            //    threadCalled[i][2] -= threadCalled[i+1][2];
            //}
            //for(int i = 0; i < superFWdepth; ++i){
            //    cout << "Depth : " << i << " took " << timeConsumePerDepth[i] << ", Called thread : "<< threadCalled[i][0] << "/"<< threadCalled[i][1] << "/"<< threadCalled[i][2] << ", sum : " << (threadCalled[i][0] + threadCalled[i][1] + threadCalled[i][2]) << ", time per thread : " << timeConsumePerDepth[i]/(threadCalled[i][0] + threadCalled[i][1] + threadCalled[i][2]) << endl;
            //}
            speedUpResultsSuperFWBFS.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsSuperFWBFS.begin(),speedUpResultsSuperFWBFS.end());
        }

        //SuperFWMET
        unsigned int superFWdepth_MET = 0;
        unsigned int superFWnum_MET = 0;
        if(SuperFWFlagMET&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            std::mutex** mutexBlock = nullptr;
            mutexBlock = new std::mutex*[nVertex];
            for(int i = 0; i < nVertex; ++i){
                mutexBlock[i] = new std::mutex[nVertex];
            }
            Supernode* root = nestedDissection(distance,nVertex,SUPERNODE_ND_MET);
            stack<pair<Supernode*,unsigned int>> operationSchedule;
            PushToScheduler(root,operationSchedule); 
            superFWdepth_MET = countDepth(root);
            superFWnum_MET = countNum(root);
            //unsigned long long timeConsumePerDepth[superFWdepth] = {0,};
            //unsigned long long threadCalled[superFWdepth][3];
            //for(int i = 0; i < superFWdepth; ++i){
            //    threadCalled[i][0] = 0;
            //    threadCalled[i][1] = 0;
            //    threadCalled[i][2] = 0;
            //}
            //printBFS(root);
            unsigned long long toPre = PAPI_get_real_usec();
            while(!operationSchedule.empty()){
                //Schedule  =>  Buffer1   (Stage1)
                //Buffer1   =>  Buffer2 
                //Buffer2   =>  free      (Stage3)
                unsigned int currentDepth = operationSchedule.top().second;
                queue<std::thread> threadStack;
                while(true){//Stage1
                    threadStack.emplace(updateEdgeSuper, operationSchedule.top().first,mutexBlock, superFWthreshold_MET,apsp9);
                    operationSchedule.pop();
                    if(operationSchedule.empty())
                        break;
                    if(operationSchedule.top().second != currentDepth)
                        break;
                }
                while(!threadStack.empty()){
                    threadStack.front().join();
                    threadStack.pop();
                }
                //timeConsumePerDepth[currentDepth] = PAPI_get_real_usec() - dTS;
                //threadCalled[currentDepth][0] = globalThreadCall[0];
                //threadCalled[currentDepth][1] = globalThreadCall[1];
                //threadCalled[currentDepth][2] = globalThreadCall[2];
            }
            //printBFS(root);
            // OutPutMatrixFile("./MatrixImage_MET.ppm",nVertex,distance,root);
            delete root;
            for(int i = 0; i < nVertex; ++i){
                    delete[] mutexBlock[i];
                }
            delete[] mutexBlock;
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "SuperFW threading MET (Post :  " << double(toPre - from)/(to - from) * 100 << "%) : " << timeConsume << "usec" << endl;
            
            // for(int i = 0; i < superFWdepth-1; ++i){
            //    threadCalled[i][0] -= threadCalled[i+1][0];
            //    threadCalled[i][1] -= threadCalled[i+1][1];
            //    threadCalled[i][2] -= threadCalled[i+1][2];
            //}
            //for(int i = 0; i < superFWdepth; ++i){
            //    cout << "Depth : " << i << " took " << timeConsumePerDepth[i] << ", Called thread : "<< threadCalled[i][0] << "/"<< threadCalled[i][1] << "/"<< threadCalled[i][2] << ", sum : " << (threadCalled[i][0] + threadCalled[i][1] + threadCalled[i][2]) << ", time per thread : " << timeConsumePerDepth[i]/(threadCalled[i][0] + threadCalled[i][1] + threadCalled[i][2]) << endl;
            //}
            speedUpResultsSuperFWMET.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsSuperFWMET.begin(),speedUpResultsSuperFWMET.end());
        }

        //Path doubling
        if(PathDubFlag&option)
        {
            unsigned long long from = PAPI_get_real_usec();
            pathDouble(nVertex, apsp8);
            unsigned long long to = PAPI_get_real_usec();
            timeConsume = (to-from);
            cout << "Path doubling : " << timeConsume << "usec" << endl;
            speedUpResultsPathDoubling.emplace_back((int)(classicTime/(double)(timeConsume) * 100));
            sort(speedUpResultsPathDoubling.begin(),speedUpResultsPathDoubling.end());
        }

        cout << "Property : |V| = " << nVertex << ", |E| = " <<nEdges ;
        if((BlockedFlag|MulBlocFlag) & option)
            cout<< ", |S| = " << sBlock;
        if(PerThreFlag & option)
            cout<< ", |T| = " << nThreads;
        if(SuperFWFlagBFS & option){
            cout<< ", SuperFWThreshold(BFS) = " << superFWthreshold_BFS << ", Depth of supernode(BFS) : " << superFWdepth_BFS << ", Num of supernode(BFS) : " << superFWnum_BFS  ;
            if(!(superFWdepth_MET & option)){
                cout << " with ";
            }
        }
        if(SuperFWFlagMET & option){
            cout<< ", SuperFWThreshold(MET) = " <<superFWdepth_MET << ", Depth of supernode(MET) : " << superFWdepth_MET << ", Num of supernode(MET) : " <<superFWnum_MET << " with " ;
        }
        cout << " " << graphName << endl;
        if(BlockedFlag&option){
            cout << "Blocked speedup : " <<
            "min : " << speedUpResultsBlocked.at(0)  << "% " <<
            "med : " << speedUpResultsBlocked.at(speedUpResultsBlocked.size() / 2)  << "% " <<
            "max : " << speedUpResultsBlocked.at(speedUpResultsBlocked.size() - 1)  << "%" << endl;
        }

        if(PerLineFlag&option){
            cout << "PerLine speedup : " <<
            "min : " << speedUpResultsPerLine.at(0)  << "% " <<
            "med : " << speedUpResultsPerLine.at(speedUpResultsPerLine.size() / 2)  << "% " <<
            "max : " << speedUpResultsPerLine.at(speedUpResultsPerLine.size() - 1)  << "%" << endl;
        }

        if(PerElemFlag&option){
            cout << "PerEle speedup : " <<
            "min : " << speedUpResultsPerElem.at(0)  << "% " <<
            "med : " << speedUpResultsPerElem.at(speedUpResultsPerElem.size() / 2)  << "% " <<
            "max : " << speedUpResultsPerElem.at(speedUpResultsPerElem.size() - 1)  << "%" << endl;
        }

        if(PerThreFlag&option){
            cout << "PerThread speedup : " <<
            "min : " << speedUpResultsPerThread.at(0)  << "% " <<
            "med : " << speedUpResultsPerThread.at(speedUpResultsPerThread.size() / 2)  << "% " <<
            "max : " << speedUpResultsPerThread.at(speedUpResultsPerThread.size() - 1)  << "%" << endl;
        }

        if((!skipClassic) && (MulBlocFlag&option)){
            cout << "BlockedThreading speedup : " <<
            "min : " << speedUpResultsBlockedqMul.at(0)  << "% " <<
            "med : " << speedUpResultsBlockedqMul.at(speedUpResultsBlockedqMul.size() / 2)  << "% " <<
            "max : " << speedUpResultsBlockedqMul.at(speedUpResultsBlockedqMul.size() - 1)  << "%" << endl;
        }
    
        if(SuperFWFlagBFS&option){
            cout << "SuperFW BFS speedup : " <<
            "min : " << speedUpResultsSuperFWBFS.at(0)  << "% " <<
            "med : " << speedUpResultsSuperFWBFS.at(speedUpResultsSuperFWBFS.size() / 2)  << "% " <<
            "max : " << speedUpResultsSuperFWBFS.at(speedUpResultsSuperFWBFS.size() - 1)  << "%" << endl;
        }

        if(SuperFWFlagMET&option){
            cout << "SuperFW MET speedup : " <<
            "min : " << speedUpResultsSuperFWMET.at(0)  << "% " <<
            "med : " << speedUpResultsSuperFWMET.at(speedUpResultsSuperFWMET.size() / 2)  << "% " <<
            "max : " << speedUpResultsSuperFWMET.at(speedUpResultsSuperFWMET.size() - 1)  << "%" << endl;
        }

        if(PathDubFlag&option){
            cout << "Path-doubling speedup : " <<
            "min : " << speedUpResultsPathDoubling.at(0)  << "% " <<
            "med : " << speedUpResultsPathDoubling.at(speedUpResultsPathDoubling.size() / 2)  << "% " <<
            "max : " << speedUpResultsPathDoubling.at(speedUpResultsPathDoubling.size() - 1)  << "%" << endl;
        }

        // Correctness handler//        
        for(int i = 0; i < nVertex; ++i){
            for(int j = 0; j < nVertex; ++j){
                bool wrongResult = false;
                string ErrorFrom;
                unsigned int criteria;
                if(skipClassic){
                    criteria = apsp6[i][j];
                }else{
                    criteria = apsp1[i][j];
                }
                if(BlockedFlag&option){
                    if(criteria != apsp2[i][j]){
                        wrongResult = true;
                        ErrorFrom = "Blocked";
                    }
                }
                if(PerLineFlag&option){
                    if(criteria != apsp3[i][j]){
                        wrongResult = true;
                        ErrorFrom = "PerLine";
                    }
                }
                if(PerElemFlag&option){
                    if(criteria != apsp4[i][j]){
                        wrongResult = true;
                        ErrorFrom = "Per element";
                    }
                }
                if(PerThreFlag&option){
                    if(criteria != apsp5[i][j]){
                        wrongResult = true;
                        ErrorFrom = "Per thread";
                    }
                }
                if((!skipClassic) && (MulBlocFlag&option)){
                    if(criteria != apsp6[i][j]){
                        wrongResult = true;
                        ErrorFrom = "Threaded block";
                    }
                }
                if(SuperFWFlagBFS&option){
                    if(criteria != apsp7[i][j]){
                        wrongResult = true;
                        ErrorFrom = "SuperFWBFS";
                    }
                }
                if(SuperFWFlagMET&option){
                    if(criteria != apsp9[i][j]){
                        wrongResult = true;
                        ErrorFrom = "SuperFWMET";
                    }
                }
                if(PathDubFlag&option){
                    if(criteria != apsp8[decomposeGPU(nVertex,i,j)]){
                        wrongResult = true;
                        ErrorFrom = "Path doubling";
                    }
                }
                if(wrongResult){
                    cout << "different result! : " << ErrorFrom<< endl;
                    exit(0);
                    
                }
            }
        }

        if(!skipClassic)
            freeResult(nVertex,apsp1);
        if(BlockedFlag&option)
            freeResult(nVertex,apsp2);
        if(PerLineFlag&option)
            freeResult(nVertex,apsp3);
        if(PerElemFlag&option)
            freeResult(nVertex,apsp4);
        if(PerThreFlag&option)
            freeResult(nVertex,apsp5);
        if(MulBlocFlag&option)
            freeResult(nVertex,apsp6);
        if(SuperFWFlagBFS&option)
            freeResult(nVertex,apsp7);
        if(PathDubFlag&option)
            freeResultGPU(apsp8);
        if(SuperFWFlagMET&option)
            freeResult(nVertex,apsp9);
        if(trial == 0){
            cout << "q to quit, r to retry, f to five trial, t to ten trial run and h to hundread trial run: ";
            cin >> terminate;
            if(terminate == 'f'){
                trial = 5;
                terminate = 'r';
            }
            if(terminate == 't'){
                trial = 10;
                terminate = 'r';
            }
            if(terminate == 'h'){
                trial = 100;
                terminate = 'r';
            }
        }
        else{
            trial--;
            terminate = 'r';
        }
    }while(terminate != 'q');
    if(distance != nullptr){   
        for(int i = 0; i < nVertex; ++i){
            delete[] distance[i];
        }
        delete[] distance;
    }
}