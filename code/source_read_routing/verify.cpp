// Independent replay: this file imports/includes no producer code.
// It consumes explicit events, not a seed or compressed routing recipe.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>
using U=uint64_t; using I=int64_t;
const U absent=UINT64_MAX;
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class T> std::vector<T> take(std::ifstream& f,U n) {
    std::vector<T> a(n); f.read(reinterpret_cast<char*>(a.data()),n*sizeof(T));
    require(bool(f),"input truncation"); return a;
}
struct Register { I value; U writer, origin; uint32_t owner; bool protectedValue; };
struct Link { uint32_t target, local, remote; };
struct Check {
    std::vector<Register> memory;
    U eid=0, scalarReads=0, scalarWrites=0, maximum=0;
    std::array<U,11> counts{};
    // Schedule fixes only the interface. Every value and writer is replayed.
    U event(U kind,U a,U b,U out,U owner,I preparation=0,U origin=absent) {
        std::array<U,8> row;
        require(std::fread(row.data(),sizeof(U),8,stdin)==8,"missing primitive event");
        require(row[0]==kind && row[1]==a && row[2]==b && row[3]==out && row[4]==owner,"event interface/schedule");
        for(int j=0;j<2;++j) {
            U key=row[1+j];
            if(key==absent) require(row[5+j]==absent,"extra read writer");
            else {
                require(key<memory.size(),"absent register read");
                require(row[5+j]==memory[key].writer && memory[key].writer<eid,"stale/forged writer");
                require(kind==5 || memory[key].owner==owner,"nonlocal read");
                ++scalarReads;
            }
        }
        I value=0; U label=origin;
        switch(kind) {
            case 0: value=0; break;
            case 1: case 2: case 10: value=preparation; break;
            case 3: require(memory[a].protectedValue,"export of mutable payload");
                    value=memory[a].value; label=memory[a].origin; break;
            case 4: require(memory[a].protectedValue && memory[a].value==0 && memory[a].origin==absent,"reset consumes payload"); value=0; break;
            case 5: {
                require(memory[b].value==0 && memory[b].origin==absent,"receiver not freshly zeroed");
                __int128 sum=__int128(memory[a].value)+memory[b].value;
                require(sum%2==0 && sum<=INT64_MAX && sum>=INT64_MIN,"mean exactness/overflow");
                value=I(sum/2); label=memory[a].origin; break;
            }
            case 6: require(memory[a].value<=INT64_MAX/2 && memory[a].value>=INT64_MIN/2,"capture overflow");
                    value=2*memory[a].value; label=memory[a].origin; break;
            case 7: require(memory[a].protectedValue && memory[a].value==0 && memory[a].origin==absent,"accumulator reset dependency"); value=2; break;
            case 8: {
                require(memory[b].protectedValue,"mutable logical read");
                require(memory[a].origin==origin,"accumulator lineage");
                __int128 sum=__int128(memory[a].value)+memory[b].value;
                require(sum<=INT64_MAX && sum>=INT64_MIN,"addition overflow"); value=I(sum); break;
            }
            case 9: value=memory[a].value+preparation; break;
            default: throw std::runtime_error("unknown opcode");
        }
        require(I(row[7])==value,"incorrect scalar operation");
        bool immutable=kind==0 || kind==2 || kind==6 || kind==9 || kind==10;
        auto put=[&](U key,U who) {
            Register next{value,eid,label,uint32_t(who),immutable};
            if(key==memory.size()) memory.push_back(next);
            else { require(key<memory.size() && !memory[key].protectedValue,"protected/absent overwrite"); memory[key]=next; }
            ++scalarWrites;
        };
        if(kind==5) { U first=memory[a].owner,second=memory[b].owner; put(a,first); put(b,second); }
        else put(out,owner);
        maximum=std::max(maximum,U(value<0?-value:value));
        ++counts[kind]; ++eid; return out;
    }
    U fresh(U kind,U owner,I prep=0,U input=absent,U origin=absent) {
        return event(kind,input,absent,memory.size(),owner,prep,origin);
    }
};
int main(int argc,char** argv) try {
    require(argc==4,"usage: verify INPUT VARIANT LOGICAL_OUTPUT");
    std::string variant=argv[2];
    require(variant=="baseline" || variant=="source" || variant=="branch" || variant=="scratch","variant");
    std::ifstream input(argv[1],std::ios::binary);
    auto h=take<U>(input,8); U q=h[0],n=h[1],c=h[2],k=h[3],centre=h[4],ports=12*c;
    require(h[7]==1 && n==q*q*q && c>=n,"header");
    auto host=take<uint32_t>(input,n), pairs=take<uint32_t>(input,4*h[5]);
    auto offsets=take<U>(input,n+1); auto neighbours=take<uint32_t>(input,h[6]);
    auto coordinates=take<I>(input,6*n); require(input.peek()==EOF,"input suffix");
    std::vector<int> site(c,-1);
    for(U s=0;s<n;++s) { require(host[s]<c && site[host[s]]==-1,"host assignment"); site[host[s]]=s; }
    std::vector<std::vector<Link>> adjacency(c);
    for(U i=0;i<h[5];++i) {
        U a=pairs[4*i],p=pairs[4*i+1],b=pairs[4*i+2],r=pairs[4*i+3];
        require(a<c && b<c && p<12 && r<12,"seam domain");
        adjacency[a].push_back({uint32_t(b),uint32_t(12*a+p),uint32_t(12*b+r)});
        adjacency[b].push_back({uint32_t(a),uint32_t(12*b+r),uint32_t(12*a+p)});
    }
    for(auto& edges:adjacency) std::sort(edges.begin(),edges.end(),[](Link x,Link y){return x.target<y.target;});
    Check v;
    // Antipodes obtained independently from the canonical icosahedral graph.
    const int axisOfPort[12]={0,1,1,0,2,3,3,2,4,5,5,4};
    const bool positive[12]={true,true,false,false,true,true,false,false,true,true,false,false};
    for(U cell=0;cell<c;++cell) for(int p=0;p<12;++p) {
        I value=0;
        if(site[cell]>=0) { I z=coordinates[6*site[cell]+axisOfPort[p]]; value=2*std::max(I(0),positive[p]?z:-z); }
        if(variant=="scratch") value+=2*(p+1);
        v.event(1,absent,absent,v.memory.size(),cell,value);
    }
    for(U cell=0;cell<c;++cell) v.fresh(0,cell);
    U accumulator=v.memory.size();
    for(U s=0;s<n;++s) v.event(1,absent,absent,v.memory.size(),host[s]);
    for(U s=0;s<n;++s) for(int a=0;a<6;++a) v.fresh(10,host[s],2*coordinates[6*s+a]);
    std::vector<U> versions(n), upcoming(n);
    std::ofstream logical(argv[3],std::ios::binary);
    auto retain=[&](U reg) {
        I value=v.memory[reg].value; require(value%2==0,"logical integer"); value/=2;
        logical.write(reinterpret_cast<char*>(&value),sizeof(I)); require(bool(logical),"logical output write");
    };
    for(U s=0;s<n;++s) { versions[s]=v.fresh(2,host[s],2*(s+1+(variant=="source" && s==centre)),absent,s); retain(versions[s]); }
    U depthSum=0,maxDepth=0,hopsPerRound=0;
    for(U layer=1;layer<=k;++layer) {
        U oldHops=v.counts[5];
        for(U t=0;t<n;++t) v.event(7,ports+host[t],absent,accumulator+t,host[t],0,layer*n+t);
        for(U s=0;s<n;++s) {
            // Independent queue traversal and backward pruning, no route manifest.
            std::vector<int> distance(c,-1),parent(c,-1),order;
            std::vector<Link> incoming(c); std::queue<uint32_t> work;
            work.push(host[s]); distance[host[s]]=0;
            while(!work.empty()) {
                uint32_t u=work.front(); work.pop(); order.push_back(u);
                for(Link edge:adjacency[u]) if(distance[edge.target]<0) {
                    distance[edge.target]=distance[u]+1; parent[edge.target]=u;
                    incoming[edge.target]=edge; work.push(edge.target);
                }
            }
            require(order.size()==c,"disconnected support");
            std::vector<bool> marked(c,false); std::vector<U> payload(c,absent);
            marked[host[s]]=true; payload[host[s]]=versions[s];
            for(U j=offsets[s];j<offsets[s+1];++j) {
                U node=host.at(neighbours[j]); maxDepth=std::max(maxDepth,U(distance[node])); depthSum+=distance[node];
                while(!marked[node]) { marked[node]=true; node=parent[node]; }
            }
            for(uint32_t node:order) {
                if(node==host[s] || !marked[node]) continue;
                Link edge=incoming[node]; U a=edge.local,b=edge.remote,local=a/12,remote=b/12;
                U archive=payload[parent[node]];
                require(archive!=absent && v.memory[archive].origin==(layer-1)*n+s,"source version chain");
                v.event(3,archive,absent,a,local);
                v.event(4,ports+remote,absent,b,remote);
                v.event(5,a,b,a,absent);
                payload[node]=v.fresh(6,remote,0,b);
                v.event(4,ports+local,absent,a,local);
                v.event(4,ports+remote,absent,b,remote);
                require(v.memory[payload[node]].value==v.memory[versions[s]].value,"transport value changed");
            }
            for(U j=offsets[s];j<offsets[s+1];++j) {
                U t=neighbours[j],data=payload[host[t]],acc=accumulator+t;
                require(data!=absent && v.memory[data].origin==(layer-1)*n+s,"logical read provenance");
                v.event(8,acc,data,acc,host[t],0,layer*n+t);
            }
        }
        for(U t=0;t<n;++t) {
            I delta=variant=="branch" && layer==1 && t==centre?2:0;
            upcoming[t]=v.fresh(9,host[t],delta,accumulator+t,layer*n+t); retain(upcoming[t]);
        }
        for(U port=0;port<ports;++port) require(v.memory[port].value==0 || v.memory[port].writer<ports,"scratch cleanup");
        versions.swap(upcoming); hopsPerRound=v.counts[5]-oldHops;
        std::cerr<<"verified layer "<<layer<<"/"<<k<<", "<<v.eid<<" events\n";
    }
    require(std::fgetc(stdin)==EOF,"extra primitive event");
    U protectedCount=0; for(const Register& r:v.memory) protectedCount+=r.protectedValue;
    std::cout<<"{\"events\":"<<v.eid<<",\"hops\":"<<v.counts[5]<<",\"register_reads\":"<<v.scalarReads
             <<",\"register_writes\":"<<v.scalarWrites<<",\"registers\":"<<v.memory.size()
             <<",\"max_abs_scaled_scalar\":"<<v.maximum<<",\"max_read_depth\":"<<maxDepth
             <<",\"total_read_depth\":"<<depthSum<<",\"tree_hops_per_round\":"<<hopsPerRound
             <<",\"protected_registers\":"<<protectedCount<<",\"mutable_registers\":"<<v.memory.size()-protectedCount
             <<",\"logical_reads\":"<<v.counts[8]<<",\"semantic_edges_checked\":"<<v.scalarReads
             <<",\"cell_struct_bytes\":"<<sizeof(Register)<<"}\n";
    return 0;
} catch(const std::exception& ex) { std::cerr<<"FAIL: "<<ex.what()<<"\n"; return 1; }
