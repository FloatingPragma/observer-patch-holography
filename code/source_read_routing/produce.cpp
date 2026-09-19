// Execute the supplied local feedback compiler. No remote scalar reads occur
// in local operations. stdout is a complete ordered tape of 64-byte events.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using U = uint64_t;
using I = int64_t;
constexpr U NONE = std::numeric_limits<U>::max();
enum Op { ZERO, SCRATCH, PUBLISH, EXPORT, RESET, MEAN, CAPTURE, START, ADD, COMMIT, ADDRESS };
struct Cell { I value; U writer; uint32_t owner; bool immutable; };
struct Edge { uint32_t to, fromPort, toPort; };
struct Engine {
    std::vector<Cell> cells;
    U count=0, reads=0, writes=0, hops=0, maxValue=0, registerHeapPeak=0;
    std::array<U,11> operations{};
    // Event columns: op, input a, input b, output, owner, writer a, writer b,
    // exact signed scalar value multiplied by two. MEAN writes a and b.
    void emit(Op op, U a, U b, U out, U owner, I value, bool immutable=false) {
        std::array<U,8> e{U(op),a,b,out,owner,a==NONE?NONE:cells.at(a).writer,
                         b==NONE?NONE:cells.at(b).writer,U(value)};
        for (U key : {a,b}) if (key!=NONE) {
            ++reads;
            if (op!=MEAN && cells.at(key).owner!=owner) throw std::runtime_error("remote read");
        }
        auto write = [&](U key, uint32_t who) {
            ++writes;
            if (key==cells.size()) {
                U oldCapacity=cells.capacity();
                cells.push_back({value,count,who,immutable});
                U allocated=cells.capacity()+(oldCapacity==cells.capacity()?0:oldCapacity);
                registerHeapPeak=std::max(registerHeapPeak,allocated*sizeof(Cell));
            }
            else {
                if (key>cells.size() || cells[key].immutable) throw std::runtime_error("immutable/absent write");
                cells[key]={value,count,who,immutable};
            }
        };
        if (op==MEAN) { write(a,cells.at(a).owner); write(b,cells.at(b).owner); }
        else write(out,uint32_t(owner));
        if (std::fwrite(e.data(),sizeof(U),8,stdout)!=8) throw std::runtime_error("event write failure");
        maxValue=std::max(maxValue,U(value<0?-value:value));
        ++count; ++operations[op];
    }
    U fresh(Op op, U owner, I value, U a=NONE) {
        U key=cells.size(); emit(op,a,NONE,key,owner,value,true); return key;
    }
    U hop(U input, U u, U v, U ports) {
        U cu=u/12, cv=v/12;
        emit(EXPORT,input,NONE,u,cu,cells.at(input).value);
        emit(RESET,ports+cv,NONE,v,cv,0);
        I sum=cells.at(u).value+cells.at(v).value;
        if (sum%2) throw std::runtime_error("half-unit representation exhausted");
        emit(MEAN,u,v,u,NONE,sum/2);
        U received=fresh(CAPTURE,cv,2*cells.at(v).value,v);
        emit(RESET,ports+cu,NONE,u,cu,0);
        emit(RESET,ports+cv,NONE,v,cv,0);
        ++hops; return received;
    }
};
template<class T> void read(std::ifstream& f, std::vector<T>& v) {
    f.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(T));
    if (!f) throw std::runtime_error("truncated input");
}
int main(int argc, char** argv) try {
    if (argc!=3) throw std::runtime_error("usage: produce INPUT baseline|source|branch|scratch");
    std::string variant=argv[2];
    if (variant!="baseline" && variant!="source" && variant!="branch" && variant!="scratch") throw std::runtime_error("variant");
    std::ifstream f(argv[1],std::ios::binary);
    std::vector<U> h(8); read(f,h);
    U q=h[0], n=h[1], c=h[2], rounds=h[3], centre=h[4], ports=12*c;
    if (h[7]!=1 || n!=q*q*q || c<n) throw std::runtime_error("input header");
    std::vector<uint32_t> host(n), pairs(4*h[5]), neighbours(h[6]);
    std::vector<U> offsets(n+1); std::vector<I> address(6*n);
    read(f,host); read(f,pairs); read(f,offsets); read(f,neighbours); read(f,address);
    if (f.peek()!=EOF) throw std::runtime_error("input suffix");
    std::vector<std::vector<Edge>> graph(c);
    for (U i=0;i<h[5];++i) {
        U a=pairs[4*i],p=pairs[4*i+1],b=pairs[4*i+2],r=pairs[4*i+3];
        graph[a].push_back({uint32_t(b),uint32_t(12*a+p),uint32_t(12*b+r)});
        graph[b].push_back({uint32_t(a),uint32_t(12*b+r),uint32_t(12*a+p)});
    }
    for (auto& row:graph) std::sort(row.begin(),row.end(),[](Edge a,Edge b){return a.to<b.to;});
    std::vector<int> siteAt(c,-1);
    for (U s=0;s<n;++s) { if (siteAt.at(host[s])!=-1) throw std::runtime_error("host alias"); siteAt[host[s]]=s; }
    Engine e;
    // Correct positive/negative antipodal placement of the six prepared axes.
    const int positive[6]={0,1,4,5,8,9}, negative[6]={3,2,7,6,11,10};
    for (U cell=0;cell<c;++cell) {
        std::array<I,12> load{};
        int s=siteAt[cell];
        if (s>=0) for (int axis=0;axis<6;++axis) {
            I z=address[6*s+axis]; load[positive[axis]]=2*std::max(I(0),z);
            load[negative[axis]]=2*std::max(I(0),-z);
        }
        for (int p=0;p<12;++p) e.emit(SCRATCH,NONE,NONE,e.cells.size(),cell,
                                      load[p]+(variant=="scratch"?2*(p+1):0));
    }
    for (U cell=0;cell<c;++cell) e.fresh(ZERO,cell,0);
    U accumulators=e.cells.size();
    for (U s=0;s<n;++s) e.emit(SCRATCH,NONE,NONE,e.cells.size(),host[s],0);
    for (U s=0;s<n;++s) for (int a=0;a<6;++a) e.fresh(ADDRESS,host[s],2*address[6*s+a]);
    std::vector<U> version(n), next(n);
    for (U s=0;s<n;++s) version[s]=e.fresh(PUBLISH,host[s],2*(s+1+(variant=="source" && s==centre)));
    std::vector<int> parent(c), depth(c), queue(c);
    std::vector<Edge> entering(c);
    std::vector<bool> used(c);
    std::vector<U> received(c);
    U maximumDepth=0, totalReadDepth=0, treeHopsPerRound=0;
    for (U layer=1;layer<=rounds;++layer) {
        U layerHops=e.hops;
        for (U s=0;s<n;++s) e.emit(START,ports+host[s],NONE,accumulators+s,host[s],2);
        for (U s=0;s<n;++s) {
            std::fill(parent.begin(),parent.end(),-1); std::fill(used.begin(),used.end(),false);
            U root=host[s], head=0, tail=1; queue[0]=root; parent[root]=root; depth[root]=0;
            while (head<tail) {
                int u=queue[head++];
                for (Edge edge:graph[u]) if (parent[edge.to]<0) {
                    parent[edge.to]=u; entering[edge.to]=edge; depth[edge.to]=depth[u]+1;
                    queue[tail++]=edge.to;
                }
            }
            if (tail!=c) throw std::runtime_error("disconnected support");
            used[root]=true; received[root]=version[s];
            for (U j=offsets[s];j<offsets[s+1];++j) {
                int u=host[neighbours[j]];
                maximumDepth=std::max(maximumDepth,U(depth[u])); totalReadDepth+=depth[u];
                while (!used[u]) { used[u]=true; u=parent[u]; }
            }
            for (U j=1;j<tail;++j) {
                int u=queue[j]; if (!used[u]) continue;
                Edge edge=entering[u];
                received[u]=e.hop(received[parent[u]],edge.fromPort,edge.toPort,ports);
            }
            for (U j=offsets[s];j<offsets[s+1];++j) {
                U target=neighbours[j], acc=accumulators+target, payload=received[host[target]];
                __int128 sum=__int128(e.cells[acc].value)+e.cells[payload].value;
                if (sum>std::numeric_limits<I>::max()) throw std::runtime_error("scalar overflow");
                e.emit(ADD,acc,payload,acc,host[target],I(sum));
            }
        }
        for (U s=0;s<n;++s) {
            U acc=accumulators+s;
            I delta=(variant=="branch" && layer==1 && s==centre)?2:0;
            next[s]=e.fresh(COMMIT,host[s],e.cells[acc].value+delta,acc);
        }
        version.swap(next); treeHopsPerRound=e.hops-layerHops;
        std::cerr<<"layer "<<layer<<"/"<<rounds<<": "<<e.count<<" events, "<<e.hops<<" hops\n";
    }
    if (std::fflush(stdout)!=0) throw std::runtime_error("event flush failure");
    U auxiliary=h.capacity()*sizeof(U)+host.capacity()*sizeof(uint32_t)
        +pairs.capacity()*sizeof(uint32_t)+neighbours.capacity()*sizeof(uint32_t)
        +offsets.capacity()*sizeof(U)+address.capacity()*sizeof(I)
        +graph.capacity()*sizeof(std::vector<Edge>)+siteAt.capacity()*sizeof(int)
        +(version.capacity()+next.capacity()+received.capacity())*sizeof(U)
        +(parent.capacity()+depth.capacity()+queue.capacity())*sizeof(int)
        +entering.capacity()*sizeof(Edge)+(used.capacity()+7)/8;
    for(const auto& row:graph) auxiliary+=row.capacity()*sizeof(Edge);
    std::cerr<<"{\"events\":"<<e.count<<",\"hops\":"<<e.hops<<",\"register_reads\":"<<e.reads
             <<",\"register_writes\":"<<e.writes<<",\"registers\":"<<e.cells.size()
             <<",\"max_abs_scaled_scalar\":"<<e.maxValue<<",\"max_read_depth\":"<<maximumDepth
             <<",\"total_read_depth\":"<<totalReadDepth<<",\"tree_hops_per_round\":"<<treeHopsPerRound
             <<",\"producer_cell_struct_bytes\":"<<sizeof(Cell)
             <<",\"register_heap_peak_requested_bytes\":"<<e.registerHeapPeak
             <<",\"auxiliary_heap_requested_bytes\":"<<auxiliary
             <<",\"wide_arithmetic_bytes\":"<<sizeof(__int128)<<",\"event_row_bytes\":"<<8*sizeof(U)<<"}\n";
    return 0;
} catch (const std::exception& ex) { std::cerr<<"FAIL: "<<ex.what()<<"\n"; return 1; }
