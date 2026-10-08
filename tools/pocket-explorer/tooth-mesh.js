'use strict';
// Parametric adult tooth families. Educational anatomy, not a measured specimen.
window.ToothMesh={build(api){
 const {a,index,upper,quad,setMaterial,center,showRoots,thirdMode}=api,c=center(a),side=Math.sign(a)||1;
 const inc=index<2,canine=index===2,pre=index===3||index===4;
 const widths=upper?[.345,.290,.320,.315,.335,.445,.415,.39]:[.275,.290,.320,.325,.345,.460,.435,.40];
 const w=widths[index],d=inc?(upper?.225:.185):canine?.265:pre?.340:(upper?.405:.385);
 const H=inc?(upper?.86:.79):canine?.94:pre?.73:.67;
 const len=Math.hypot(2.55*Math.cos(a),3.2*Math.sin(a)),tx=2.55*Math.cos(a)/len,tz=-3.2*Math.sin(a)/len,nx=-tz,nz=tx;
 function P(x,y,z){if(index===7){if(thirdMode==='partial')y-=.30;if(thirdMode==='impacted'){let old=x;x=.57*x+.82*y;y=-.82*old+.57*y-.60;}}return[c[0]+tx*x+nx*z,y+.025*Math.abs(a),c[2]+tz*x+nz*z];}
 const sp=(x,n)=>Math.sign(x)*Math.abs(x)**n;
 function outline(t){const x=sp(Math.cos(t),inc?.68:canine?.84:pre?.87:.61),z=sp(Math.sin(t),inc?.78:canine?.92:pre?.90:.67);const skew=upper&&!inc&&!canine&&!pre?.10*z:0;return[w*(x+skew)*(1-.045*side*x),d*z*(1+.035*x)];}
 function top(x,z){const u=x/w,v=z/d;if(inc)return H-.046*u*u-.014*Math.cos(u*5.5);if(canine)return .69+.24*Math.exp(-3.3*u*u-3*v*v)-.035*Math.max(0,-v);
 const cusps=pre?[[0,.54,.20],[0,-.53,upper?.17:.12]]:upper?[[-.47,.50,.16],[.48,.48,.135],[-.46,-.48,.19],[.51,-.47,.14]]:[[-.53,.51,.17],[.43,.52,.16],[-.50,-.52,.20],[.45,-.54,.18],...(index===5?[[.85,.04,.075]]:[])];
 let h=pre?.52:.46;for(const [cx,cz,height]of cusps)h+=height*Math.exp(-4.5*((u-cx)**2+(v-cz)**2));
 h-=.027*Math.exp(-100*v*v)*Math.exp(-1.8*u*u);if(!pre)h-=.019*Math.exp(-90*(u+.07*Math.sin(v*4))**2)*Math.exp(-2*v*v);
 // Marginal ridges bound the occlusal table without separate spherical cusps.
 h+=.024*Math.exp(-60*(Math.abs(u)-.83)**2)*(1-.4*v*v);return h;
 }
 function shape(r,v){const f=.74+.30*Math.sin(v*Math.PI*.77),fz=inc?(.94-.71*v*v):canine?(.78+.24*Math.sin(v*Math.PI*.85)):f;let x=r[0]*f,z=r[1]*fz;
 if(inc){z+=.035*v;z+=.008*Math.cos(x/w*8)*Math.max(0,r[1]/d)*Math.sin(v*Math.PI);if(r[1]<0)z+=.024*Math.sin(v*Math.PI);}
 const y=-.09+v*(top(x,z)+.09);return [x,y,z];}
 setMaterial('enamel');const n=56,rows=18,rings=18;
 for(let j=0;j<n;j++){const b=outline(j/n*Math.PI*2),e=outline((j+1)/n*Math.PI*2);for(let k=0;k<rows;k++){let v=k/rows,vn=(k+1)/rows;const col=[230+12*v,216+22*v,188+37*v];quad(P(...shape(b,v)),P(...shape(b,vn)),P(...shape(e,vn)),P(...shape(e,v)),col);}
 const be=shape(b,1),ee=shape(e,1);for(let k=0;k<rings;k++){let r=k/rings,rn=(k+1)/rings;const Q=(p,f)=>P(p[0]*f,top(p[0]*f,p[2]*f),p[2]*f);const shade=.987-.02*(1-r);quad(Q(be,r),Q(ee,r),Q(ee,rn),Q(be,rn),[244*shade,240*shade,222*shade]);}}
 if(!showRoots)return;
 setMaterial('root');const branches=inc||canine?[[0,0]]:pre?(upper&&index===3?[[-.06,.10],[.06,-.10]]:[[0,0]]):upper?[[-.18,.15],[.17,.17],[0,-.21]]:[[-.21,0],[.21,0]];
 // Cervical trunk bridges crown and furcation instead of isolated cone roots.
 const trunk=branches.length>1?.30:.05;
 for(let k=0;k<6;k++)for(let j=0;j<32;j++){const Q=(v,t)=>P(w*.74*(1-.22*v)*Math.cos(t),-.09-v*trunk,d*.80*(1-.22*v)*Math.sin(t));quad(Q(k/6,j/32*6.283185),Q(k/6,(j+1)/32*6.283185),Q((k+1)/6,(j+1)/32*6.283185),Q((k+1)/6,j/32*6.283185),[208,190,153]);}
 for(const [rx,rz]of branches)for(let k=0;k<24;k++)for(let j=0;j<24;j++){const Q=(v,t)=>{const radius=(branches.length===1?(inc?.15:.18):.14)*Math.pow(Math.max(0,1-v),.64)+.003;return P(rx*(.65+v*.75)+radius*Math.cos(t)+side*.04*v*v,-.09-trunk-v*(canine?1.38:1.08),rz*(.60+v*.70)+radius*.80*Math.sin(t));};quad(Q(k/24,j/24*6.283185),Q(k/24,(j+1)/24*6.283185),Q((k+1)/24,(j+1)/24*6.283185),Q((k+1)/24,j/24*6.283185),[205,186,148]);}
}};
