#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using U=uint64_t;
struct RNG {
 U x; explicit RNG(U seed):x(seed){}
 U next(){ U z=(x+=0x9e3779b97f4a7c15ULL); z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
 z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31); }
 U below(U n){if(!n)throw std::runtime_error("zero range"); U t=-n%n,z;do{z=next();}while(z<t);return z%n;}
};
U add_checked(U a,U b){if(a>std::numeric_limits<U>::max()-b)throw std::runtime_error("integer overflow");return a+b;}
U mul_checked(U a,U b){if(b && a>std::numeric_limits<U>::max()/b)throw std::runtime_error("integer overflow");return a*b;}
struct Node {int a,b; U chains,birth_tau,birth_step,write_tau=0,write_step=0; std::vector<U> anc;};
struct State {
 std::vector<Node> nodes; std::unordered_map<U,uint32_t> pairs;
 U F=0,P=0,tau=0,t=0,work=0,forced=0,relieved=0; int H;
 explicit State(int h):H(h){nodes.push_back({-1,-1,1,0,0,0,0,{1}});nodes.push_back({-1,-1,1,0,0,0,0,{2}});}
 U weight(size_t i)const{return mul_checked(2,nodes[i].chains)-2;}
 void step(const std::vector<U>& drawn,std::vector<uint32_t>& W,U &C,U &v){
  auto n=nodes.size(); U q=std::min<U>(3,F+n); if(drawn.size()!=q)throw std::runtime_error("bad q");
  W.clear(); for(U z:drawn){if(z>=F+n)throw std::runtime_error("bad draw"); if(z<n)W.push_back(z);}
  auto sd=drawn;std::sort(sd.begin(),sd.end());if(std::adjacent_find(sd.begin(),sd.end())!=sd.end())throw std::runtime_error("duplicate draw");
  std::sort(W.begin(),W.end());U f=q-W.size(),fm=F-f,pm=add_checked(P,2*f);
  U quota=2*((std::max<U>(1,pm/6)+1)/2); v=(fm>=3 && pm>=6)?std::min({quota,U(H),fm,pm}):0;
  std::vector<std::array<uint32_t,2>> batch;
  for(size_t i=0;i<W.size();++i)for(size_t j=i+1;j<W.size();++j){U key=(U(W[i])<<32)|W[j];if(!pairs.count(key))batch.push_back({W[i],W[j]});}
  // Sequential first/repeat charging equals complete-burst formula; all parents preexist.
  C=0; U newtau=add_checked(tau,W.size()),newt=t+1;
  for(auto uv:batch){uint32_t a=uv[0],b=uv[1],z=nodes.size();auto &va=nodes[a].anc;auto &vb=nodes[b].anc;
   size_t nw=z/64+1;std::vector<U> an(nw,0);
   for(size_t k=0;k<va.size();++k){an[k]|=va[k];} for(size_t k=0;k<vb.size();++k){an[k]|=vb[k];}
   U cost=0;
   for(size_t k=0;k<nw;++k){U mask=an[k];while(mask){int bit=__builtin_ctzll(mask);size_t w=k*64+bit;mask&=mask-1;
    cost=add_checked(cost,mul_checked(nodes[w].write_step?2:11,weight(w)));}}
   C=add_checked(C,cost);
   for(uint32_t p:uv)if(!nodes[p].write_step){nodes[p].write_step=newt;nodes[p].write_tau=newtau;}
   U chains=add_checked(nodes[a].chains,nodes[b].chains);an[z/64]|=U(1)<<(z%64);
   nodes.push_back({int(a),int(b),chains,newtau,newt,0,0,std::move(an)});pairs[(U(a)<<32)|b]=z;
  }
  F=add_checked(fm-v,C);P=pm-v;t=newt;tau=newtau;
  work=add_checked(work,C);forced=add_checked(forced,f);relieved=add_checked(relieved,v);
  if(F!=work-forced-relieved || P!=2*forced-relieved || tau+forced!=3*t-1)throw std::runtime_error("budget");
 }
};
// Event records: six uint64 (step,tau,n_before,F_before,P_before,cost),
// four uint32 (served count,served0/1/2), for 64 bytes. Little-endian hosts only.
int main_standard(int argc,char**argv){try{
 if(argc<6)throw std::runtime_error("usage: native35 prefix H seed target_tau max_steps [audit_steps]");
 std::string prefix=argv[1];int H=std::stoi(argv[2]);U seed=std::stoull(argv[3]),target=std::stoull(argv[4]),maxsteps=std::stoull(argv[5]);
 U audit=argc>6?std::stoull(argv[6]):0; if(H!=0&&H!=6)throw std::runtime_error("H must 0 or6");
 unsigned endian=1;if(*(char*)&endian!=1)throw std::runtime_error("little endian required");
 std::ofstream events(prefix+".events.bin",std::ios::binary),out(prefix+".nodes.tsv"),au(prefix+".audit.tsv");
 if(!events||!out||!au)throw std::runtime_error("output paths");
 State st(H);RNG rng(seed);std::vector<uint32_t> W;std::vector<U> draw;U lastTau=0;
 auto start=std::chrono::steady_clock::now();
 while(st.tau<target && st.t<maxsteps && st.nodes.size()<50000){
  U n=st.nodes.size(),F=st.F,P=st.P,R=n+F,q=std::min<U>(3,R);draw.clear();
  for(U j=0;j<q;++j){U x;do{x=rng.below(R);}while(std::find(draw.begin(),draw.end(),x)!=draw.end());draw.push_back(x);}
  U C,v;st.step(draw,W,C,v);
  if(!W.empty()){
   U fields[6]={st.t,st.tau,n,F,P,C};uint32_t w[4]={uint32_t(W.size()),UINT32_MAX,UINT32_MAX,UINT32_MAX};
   for(size_t k=0;k<W.size();++k){w[k+1]=W[k];} events.write((char*)fields,sizeof(fields));events.write((char*)w,sizeof(w));
  }
  if(st.t<=audit){au<<st.t<<'\t'<<n<<'\t'<<F<<'\t'<<P;for(U x:draw)au<<'\t'<<x;for(U j=q;j<3;++j)au<<"\t-1";
   au<<'\t'<<st.nodes.size()<<'\t'<<st.F<<'\t'<<st.P<<'\t'<<st.tau<<'\t'<<C<<'\t'<<v<<'\n';}
  if(st.tau/65536!=lastTau/65536){std::cerr<<prefix<<" tau="<<st.tau<<" n="<<st.nodes.size()<<" steps="<<st.t<<'\n';lastTau=st.tau;}
 }
 for(size_t i=0;i<st.nodes.size();++i){auto &a=st.nodes[i];out<<i<<'\t'<<a.a<<'\t'<<a.b<<'\t'<<a.chains<<'\t'<<a.birth_tau<<'\t'<<a.birth_step<<'\t'<<a.write_tau<<'\t'<<a.write_step<<'\n';}
 std::ofstream meta(prefix+".meta.json");double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 meta<<"{\"seed\":"<<seed<<",\"H\":"<<H<<",\"target_tau\":"<<target<<",\"tau\":"<<st.tau<<",\"steps\":"<<st.t<<",\"n\":"<<st.nodes.size()<<",\"F\":"<<st.F<<",\"P\":"<<st.P<<",\"total_work\":"<<st.work<<",\"forced_served\":"<<st.forced<<",\"relief_removed\":"<<st.relieved<<",\"rng_state\":"<<rng.x<<",\"complete\":"<<(st.tau>=target?"true":"false")<<",\"runtime_seconds\":"<<secs<<"}\n";
 std::cout<<"tau="<<st.tau<<" t="<<st.t<<" n="<<st.nodes.size()<<" seconds="<<secs<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

int fork_main(int argc,char**argv){
 if(argc!=7)throw std::runtime_error("--fork nodes_header output seed count max_steps");
 std::ifstream input(argv[2]);std::string output=argv[3];U seed=std::stoull(argv[4]),reps=std::stoull(argv[5]),maxsteps=std::stoull(argv[6]);
 int H;U F,P,t,tau,n;input>>H>>F>>P>>t>>tau>>n;
 State base(H);base.nodes.clear();base.pairs.clear();base.F=F;base.P=P;base.t=t;base.tau=tau;
 base.forced=3*t-1-tau;base.relieved=2*base.forced-P;base.work=F+base.forced+base.relieved;
 for(U i=0;i<n;++i){Node a;input>>a.a>>a.b>>a.chains>>a.birth_tau>>a.birth_step>>a.write_tau>>a.write_step;
  a.anc.assign(i/64+1,0);if(i>=2){for(int p:{a.a,a.b}){for(size_t k=0;k<base.nodes[p].anc.size();++k)a.anc[k]|=base.nodes[p].anc[k];}base.pairs[(U(a.a)<<32)|a.b]=i;}
  a.anc[i/64]|=U(1)<<(i%64);base.nodes.push_back(std::move(a));
 }
 int ka,kb;input>>ka;std::vector<int>A(ka);for(auto&x:A)input>>x;input>>kb;std::vector<int>B(kb);for(auto&x:B)input>>x;if(!input)throw std::runtime_error("fork input malformed");
 std::ofstream out(output);if(!out)throw std::runtime_error("fork output");
 for(U k=0;k<reps;++k){State st=base;RNG rng(seed+k);U na=0,nb=0,un=0;bool gotg=false,gotl=false;
  auto save=[&](const std::string&mode,bool done){out<<k<<'\t'<<seed+k<<'\t'<<mode<<'\t'<<int(done)<<'\t'<<st.t-t<<'\t'<<st.tau-tau<<'\t'<<na<<'\t'<<nb<<'\t'<<un<<'\t'<<st.nodes.size()<<'\t'<<st.F<<'\n';};
  while((!gotg||!gotl)&&st.t-t<maxsteps&&st.tau-tau<262144&&st.nodes.size()<50000){
   U R=st.F+st.nodes.size(),q=std::min<U>(3,R);std::vector<U>draw;
   for(U j=0;j<q;++j){U x;do{x=rng.below(R);}while(std::find(draw.begin(),draw.end(),x)!=draw.end());draw.push_back(x);}
   U C,v;std::vector<uint32_t>W;st.step(draw,W,C,v);
   for(int z:W){bool a=std::find(A.begin(),A.end(),z)!=A.end(),b=std::find(B.begin(),B.end(),z)!=B.end();na+=a;nb+=b;un+=a||b;}
   if(!gotg&&st.tau-tau>=8192){save("global8192",true);gotg=true;}
   if(!gotl&&un>=32){save("local32",true);gotl=true;}
  }
  if(!gotg){save("global8192",false);} if(!gotl){save("local32",false);}
 }
 return 0;
}
int main(int argc,char**argv){try{if(argc>1&&std::string(argv[1])=="--fork")return fork_main(argc,argv);return main_standard(argc,argv);}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
