/* Procedural oral soft tissues. Illustrative geometry in jaw-local model units. */
(function(root){
'use strict';
const TAU=Math.PI*2, mix=(a,b,t)=>a+(b-a)*t;
function patch(api,nu,nv,point,color){for(let i=0;i<nu;i++)for(let j=0;j<nv;j++){const u=i/nu,v=j/nv;api.quad(point(u,v),point((i+1)/nu,v),point((i+1)/nu,(j+1)/nv),point(u,(j+1)/nv),typeof color==='function'?color(u,v):color);}}
function tube(api,points,r,color){const n=points.length-1;for(let i=0;i<n;i++){const c=points[i],d=points[i+1],t=d.map((x,k)=>x-c[k]),len=Math.hypot(...t)||1;t.forEach((x,k)=>t[k]=x/len);let b=Math.abs(t[1])<.9?[0,1,0]:[1,0,0];let p=[t[1]*b[2]-t[2]*b[1],t[2]*b[0]-t[0]*b[2],t[0]*b[1]-t[1]*b[0]],pl=Math.hypot(...p);p=p.map(x=>x/pl);b=[t[1]*p[2]-t[2]*p[1],t[2]*p[0]-t[0]*p[2],t[0]*p[1]-t[1]*p[0]];const q=(c,a)=>c.map((x,k)=>x+r*(p[k]*Math.cos(a)+b[k]*Math.sin(a)));for(let j=0;j<10;j++)api.quad(q(c,j/10*TAU),q(d,j/10*TAU),q(d,(j+1)/10*TAU),q(c,(j+1)/10*TAU),color);}}
function gingiva(api){api.setMaterial('gum');const angles=[.12,.36,.61,.85,1.07,1.28,1.55],col=api.upper?[186,107,115]:[178,94,106];
// Closed, broad alveolar ridge: cervical contours taper into mucosal wall, not a circular tube.
const p=(u,v)=>{const a=-1.77+3.54*u,c=api.center(a),t=v*TAU;let nearest=Math.min(...angles.map(b=>Math.abs(Math.abs(a)-b)));const papilla=.085*Math.min(1,nearest/.10),top=-.105+papilla;const lateral=Math.cos(t),roof=Math.max(0,Math.sin(t)),width=.36+.06*Math.abs(a)/1.77;return[c[0]+Math.sin(a)*width*lateral,top-.30+.30*Math.sin(t)-.04*(1-roof),c[2]+Math.cos(a)*width*lateral];};patch(api,124,24,p,(u,v)=>{const f=.94+.06*Math.sin(v*TAU);return col.map(x=>x*f);});
}
function palate(api){api.setMaterial('mucosa');
// Vault rises superiorly in reflected upper coordinates. Edges join lingual ridge.
const roof=(u,v)=>{const z=mix(-.65,2.86,v),w=2.14*Math.sqrt(Math.max(.03,1-((z-.12)/3.05)**2)),x=(u*2-1)*w,arch=Math.pow(Math.max(0,1-(x/w)**2),.7);return[x,-.18-.66*arch*Math.sin(Math.PI*(.16+.70*v)),z];};
patch(api,44,48,roof,(u,v)=>[177+8*v,105+10*v,113+10*v]);
// Low-relief rugae follow the anterior vault, kept subordinate to anatomy.
for(let side of [-1,1])for(let k=0;k<5;k++){const z0=1.23+k*.245,pts=[];for(let j=0;j<=22;j++){const s=j/22,x=side*(.08+s*(.94-k*.09)),z=z0-.19*Math.sin(s*Math.PI/2),v=(z+.65)/3.51,w=2.14*Math.sqrt(Math.max(.03,1-((z-.12)/3.05)**2)),u=(x/w+1)/2;let p=roof(u,v);p[1]+=.014;pts.push(p);}tube(api,pts,.018,[190,118,125]);}
const raphe=[];for(let j=0;j<=28;j++){let p=roof(.5,.16+j/28*.76);p[1]+=.012;raphe.push(p);}tube(api,raphe,.011,[187,116,121]);
// A continuous posterior drape overlaps the hard-palate edge.
patch(api,40,22,(u,v)=>{const x=(2*u-1)*mix(2.10,1.25,v),z=-.65-.70*v,y=mix(-.18-.318*Math.pow(Math.max(0,1-(2*u-1)**2),.7),.17,v);return[x,y,z];},[163,83,96]);
api.ellipsoid([0,.12,-1.30],[.105,.29,.13],[169,87,100]);
}
function tongue(api){api.setMaterial('tongue');
// Rounded superelliptic body, a broad root and restrained tapered tip.
patch(api,58,34,(u,v)=>{const z=-.87+3.18*u,t=v*TAU,end=Math.pow(Math.max(0,Math.sin(Math.PI*u)),.40),width=(1.60-.39*u)*end,x=width*Math.cos(t),top=Math.sin(t),groove=.047*Math.exp(-x*x/.028)*Math.max(0,top)*Math.sin(Math.PI*u);return[x,-.29+.50*end*top-.035*u-groove,z];},(u,v)=>{const d=Math.max(0,Math.sin(v*TAU));return[153+22*d,76+12*d,94+13*d];});
// Floor of mouth is a shallow crescent rather than another floating oval.
api.setMaterial('mucosa');patch(api,58,12,(u,v)=>{const a=-1.65+3.30*u,c=api.center(a),r=mix(.60,.94,v);return[c[0]*r,-.60+.23*v,c[2]*r-.08];},[139,72,87]);
}
function cheeks(api){api.setMaterial('mucosa');for(let side of [-1,1])patch(api,42,20,(u,v)=>{const a=side*(.67+1.20*u),c=api.center(a),out=.46+.23*Math.sin(Math.PI*v);return[c[0]+Math.sin(a)*out,-.67+1.25*v,c[2]+Math.cos(a)*out];},[147,73,88]);}
function lips(api){api.setMaterial('lip');
// Retracted vermilion ribbon, with a shallow central bow, left open for inspection.
patch(api,70,18,(u,v)=>{const a=-1.03+2.06*u,c=api.center(a),t=v*TAU,bow=.045*Math.cos(a*6)*Math.exp(-a*a/.35),r=.42+.14*Math.cos(t);return[c[0]*1.14+Math.sin(a)*r,-.33+bow+.15*Math.sin(t),c[2]+Math.cos(a)*r];},api.upper?[171,80,96]:[179,87,104]);}
function posterior(api){api.setMaterial('mucosa');
// Rounded concave posterior envelope. Its lateral rim returns forward behind
// the gum pads; the low rim is buried in the tongue instead of hanging free.
patch(api,64,42,(u,v)=>{const t=(u*2-1)*Math.PI/2,y=-.06+2.13*v;
const width=1.65+.22*Math.sin(Math.PI*v),x=width*Math.sin(t);
const z=-1.64+.78*Math.pow(Math.abs(Math.sin(t)),2)+.85*Math.pow(v,5);
// Two subtle integral ridges suggest the pillars without disconnected tubing.
const ridge=.065*Math.exp(-Math.pow((Math.abs(t)-1.03)/.13,2))*Math.sin(Math.PI*v);
return[x,y,z+ridge];},(u,v)=>{const rim=Math.pow(Math.abs(u*2-1),2);return[126+17*rim,58+10*rim,75+12*rim];});
// Broad rounded return over the superior rim connects the posterior drape.
patch(api,48,16,(u,v)=>{const x=(u*2-1)*1.63;
return[x,mix(-.17,-.06,v),mix(-.87,-1.64+.78*Math.pow(x/1.65,2),v)];},[150,73,88]);
}
root.MouthGeometry={gingiva,palate,tongue,cheeks,lips,posterior};
})(typeof window!=='undefined'?window:globalThis);
