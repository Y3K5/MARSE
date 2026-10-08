import {ANATOMY,MODEL,type TissueId} from './anatomy';
export type Mode='molar'|'mouth';
export interface SceneState {schema:1;model:string;version:string;mode:Mode;section:boolean;cut:number;explode:number;exposure:number;lighting:'studio'|'clinical'|'rim';labels:boolean;quality:'balanced'|'cinematic';selected:TissueId;hidden:TissueId[];camera:number[];target:number[];}
export function validateScene(value:unknown):SceneState {
 if(!value||typeof value!=='object')throw Error('Scene must be a JSON object.');
 const s=value as SceneState;
 if(s.schema!==1||s.model!==MODEL.id||s.version!==MODEL.version)throw Error('Scene schema or model version does not match this viewer.');
 const ids=new Set(ANATOMY.map(a=>a.id));
 const bounded=(x:unknown,min:number,max:number)=>typeof x==='number'&&Number.isFinite(x)&&x>=min&&x<=max;
 if(!['molar','mouth'].includes(s.mode)||!['studio','clinical','rim'].includes(s.lighting)||!['balanced','cinematic'].includes(s.quality)||!ids.has(s.selected))throw Error('Invalid scene setting.');
 if(!bounded(s.cut,-5,5)||!bounded(s.explode,0,1)||!bounded(s.exposure,.3,2)||typeof s.section!=='boolean'||typeof s.labels!=='boolean')throw Error('Invalid scene controls.');
 if(!Array.isArray(s.hidden)||s.hidden.some(id=>!ids.has(id))||new Set(s.hidden).size!==s.hidden.length)throw Error('Invalid tissue visibility.');
 for(const v of [s.camera,s.target])if(!Array.isArray(v)||v.length!==3||v.some(x=>!bounded(x,-500,500)))throw Error('Invalid camera coordinates.');
 if(Math.hypot(...s.camera.map((x,i)=>x-s.target[i]))<.1)throw Error('Camera must be separated from its target.');
 return structuredClone(s);
}
