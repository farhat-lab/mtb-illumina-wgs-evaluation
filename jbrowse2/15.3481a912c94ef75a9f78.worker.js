this.webpackChunk([15],{1011:function(e,t){e.exports=function(e,t){if(!(e instanceof t))throw new TypeError("Cannot call a class as a function")}},1012:function(e,t){function r(e,t){for(var r=0;r<t.length;r++){var n=t[r];n.enumerable=n.enumerable||!1,n.configurable=!0,"value"in n&&(n.writable=!0),Object.defineProperty(e,n.key,n)}}e.exports=function(e,t,n){return t&&r(e.prototype,t),n&&r(e,n),e}},1027:function(e,t,r){"use strict";var n=r(1100),a=r(1038),i=a.unzip,u=a.unzipChunk,s=a.unzipChunkSlice;e.exports={BgzfFilehandle:n,unzip:i,unzipChunk:u,unzipChunkSlice:s}},1038:function(e,t,r){"use strict";(function(t){var n=r(997),a=n(r(998)),i=n(r(999)),u=r(683),s=(0,r(682).promisify)(u.gunzip),o=r(375),c=o.Z_SYNC_FLUSH,f=o.Inflate;function l(e){return h.apply(this,arguments)}function h(){return(h=(0,i.default)(a.default.mark((function e(r){var n,i,u,s,o,l,h;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:e.prev=0,i=0,u=0,s=[];case 4:if(l=r.slice(i),o=new f,n=o.strm,o.push(l,c),!o.err){e.next=11;break}throw new Error(o.msg);case 11:i+=n.next_in,s[u]=t.from(o.result),u+=1;case 14:if(n.avail_in){e.next=4;break}case 15:return h=t.concat(s),e.abrupt("return",h);case 19:if(e.prev=19,e.t0=e.catch(0),!e.t0.message.match(/incorrect header check/)){e.next=23;break}throw new Error("problem decompressing block: incorrect gzip header check");case 23:case"end":return e.stop()}}),e,null,[[0,19]])})))).apply(this,arguments)}function p(){return(p=(0,i.default)(a.default.mark((function e(r){var n,i,u,s,o,l,h,p,d,v;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:e.prev=0,i=0,u=0,s=[],o=[],l=[];case 6:if(h=r.slice(i),p=new f,n=p.strm,p.push(h,c),!p.err){e.next=12;break}throw new Error(p.msg);case 12:d=t.from(p.result),s.push(d),o.push(i),l.push(u),i+=n.next_in,u+=d.length;case 18:if(n.avail_in){e.next=6;break}case 19:return v=t.concat(s),e.abrupt("return",{buffer:v,cpositions:o,dpositions:l});case 23:if(e.prev=23,e.t0=e.catch(0),!e.t0.message.match(/incorrect header check/)){e.next=27;break}throw new Error("problem decompressing block: incorrect gzip header check");case 27:case"end":return e.stop()}}),e,null,[[0,23]])})))).apply(this,arguments)}function d(){return(d=(0,i.default)(a.default.mark((function e(r,n){var i,u,s,o,l,h,p,d,v,b,x,k;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:e.prev=0,u=n.minv.blockPosition,s=n.minv.dataPosition,o=[],l=[],h=[];case 6:if(p=r.slice(u-n.minv.blockPosition),d=new f,i=d.strm,d.push(p,c),!d.err){e.next=12;break}throw new Error(d.msg);case 12:if(v=t.from(d.result),o.push(v),b=v.length,l.push(u),h.push(s),1===o.length&&n.minv.dataPosition&&(o[0]=o[0].slice(n.minv.dataPosition),b=o[0].length),x=u,u+=i.next_in,s+=b,!(x>=n.maxv.blockPosition)){e.next=26;break}return o[o.length-1]=o[o.length-1].slice(0,n.maxv.blockPosition===n.minv.blockPosition?n.maxv.dataPosition-n.minv.dataPosition+1:n.maxv.dataPosition+1),l.push(u),h.push(s),e.abrupt("break",27);case 26:if(i.avail_in){e.next=6;break}case 27:return k=t.concat(o),e.abrupt("return",{buffer:k,cpositions:l,dpositions:h});case 31:if(e.prev=31,e.t0=e.catch(0),!e.t0.message.match(/incorrect header check/)){e.next=35;break}throw new Error("problem decompressing block: incorrect gzip header check");case 35:case"end":return e.stop()}}),e,null,[[0,31]])})))).apply(this,arguments)}e.exports={unzip:l,unzipChunk:function(e){return p.apply(this,arguments)},unzipChunkSlice:function(e,t){return d.apply(this,arguments)},nodeUnzip:function(e){return s(e,{finishFlush:(u.constants||u).Z_SYNC_FLUSH})},pakoUnzip:l}}).call(this,r(74).Buffer)},1039:function(e,t,r){"use strict";var n=r(997),a=n(r(998)),i=n(r(999)),u=n(r(1011)),s=n(r(1012)),o=void 0,c=function(){function e(t){(0,u.default)(this,e),this.fdPromise=o.open(t,"r"),this.path=t}return(0,s.default)(e,[{key:"read",value:function(){var e=(0,i.default)(a.default.mark((function e(t,r,n,i){var u,s;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.fdPromise;case 2:return u=e.sent,e.next=5,o.read(u,t,r,n,i);case 5:return s=e.sent,e.abrupt("return",s);case 7:case"end":return e.stop()}}),e,this)})));return function(t,r,n,a){return e.apply(this,arguments)}}()},{key:"stat",value:function(){var e=(0,i.default)(a.default.mark((function e(){var t;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.fdPromise;case 2:return t=e.sent,e.abrupt("return",o.fstat(t));case 4:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()}]),e}();e.exports=c},1100:function(e,t,r){"use strict";(function(t){var n=r(997),a=n(r(1101)),i=n(r(998)),u=n(r(999)),s=n(r(1011)),o=n(r(1012)),c=r(1038).unzip,f=r(1039),l=r(1107),h=function(){function e(t){var r=t.filehandle,n=t.path,a=t.gziFilehandle,i=t.gziPath;if((0,s.default)(this,e),r)this.filehandle=r;else{if(!n)throw new TypeError("either filehandle or path must be defined");this.filehandle=new f(n)}if(!a&&!i&&!n)throw new TypeError("either gziFilehandle or gziPath must be defined");this.gzi=new l({filehandle:a,path:a||i||!n?"".concat(n,".gzi"):i})}return(0,o.default)(e,[{key:"stat",value:function(){var e=(0,u.default)(i.default.mark((function e(){var t;return i.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.filehandle.stat();case 2:return t=e.sent,e.t0=Object,e.t1=t,e.next=7,this.getUncompressedFileSize();case 7:return e.t2=e.sent,e.t3=void 0,e.t4=void 0,e.t5={size:e.t2,blocks:e.t3,blksize:e.t4},e.abrupt("return",e.t0.assign.call(e.t0,e.t1,e.t5));case 12:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"getUncompressedFileSize",value:function(){var e=(0,u.default)(i.default.mark((function e(){var r,n,u,s,o,c,f,l;return i.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.gzi.getLastBlock();case 2:return r=e.sent,n=(0,a.default)(r,2),u=n[1],e.next=7,this.filehandle.stat();case 7:return s=e.sent,o=s.size,c=t.allocUnsafe(4),e.next=12,this.filehandle.read(c,0,4,o-28-4);case 12:if(f=e.sent,4===f.bytesRead){e.next=16;break}throw new Error("read error");case 16:return l=c.readUInt32LE(0),e.abrupt("return",u+l);case 18:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"_readAndUncompressBlock",value:function(){var e=(0,u.default)(i.default.mark((function e(t,r,n){var u,s,o,f,l,h,p;return i.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:if(u=(0,a.default)(r,1),s=u[0],o=(0,a.default)(n,1),f=o[0],l=f){e.next=7;break}return e.next=6,this.filehandle.stat();case 6:l=e.sent.size;case 7:return h=l-s,e.next=10,this.filehandle.read(t,0,h,s);case 10:return e.next=12,c(t.slice(0,h));case 12:return p=e.sent,e.abrupt("return",p);case 14:case"end":return e.stop()}}),e,this)})));return function(t,r,n){return e.apply(this,arguments)}}()},{key:"read",value:function(){var e=(0,u.default)(i.default.mark((function e(r,n,u,s){var o,c,f,l,h,p,d,v,b,x;return i.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.gzi.getRelevantBlocksForRead(u,s);case 2:o=e.sent,c=t.allocUnsafe(65536),f=n,l=0,h=0;case 7:if(!(h<o.length-1)){e.next=18;break}return e.next=10,this._readAndUncompressBlock(c,o[h],o[h+1]);case 10:p=e.sent,d=(0,a.default)(o[h],2),v=d[1],b=v>=s?0:s-v,x=Math.min(s+u,v+p.length)-v,b>=0&&b<p.length&&(p.copy(r,f,b,x),f+=x-b,l+=x-b);case 15:h+=1,e.next=7;break;case 18:return e.abrupt("return",{bytesRead:l,buffer:r});case 19:case"end":return e.stop()}}),e,this)})));return function(t,r,n,a){return e.apply(this,arguments)}}()}]),e}();e.exports=h}).call(this,r(74).Buffer)},1101:function(e,t,r){var n=r(1102),a=r(1103),i=r(1104),u=r(1106);e.exports=function(e,t){return n(e)||a(e,t)||i(e,t)||u()}},1102:function(e,t){e.exports=function(e){if(Array.isArray(e))return e}},1103:function(e,t){e.exports=function(e,t){if("undefined"!==typeof Symbol&&Symbol.iterator in Object(e)){var r=[],n=!0,a=!1,i=void 0;try{for(var u,s=e[Symbol.iterator]();!(n=(u=s.next()).done)&&(r.push(u.value),!t||r.length!==t);n=!0);}catch(o){a=!0,i=o}finally{try{n||null==s.return||s.return()}finally{if(a)throw i}}return r}}},1104:function(e,t,r){var n=r(1105);e.exports=function(e,t){if(e){if("string"===typeof e)return n(e,t);var r=Object.prototype.toString.call(e).slice(8,-1);return"Object"===r&&e.constructor&&(r=e.constructor.name),"Map"===r||"Set"===r?Array.from(e):"Arguments"===r||/^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(r)?n(e,t):void 0}}},1105:function(e,t){e.exports=function(e,t){(null==t||t>e.length)&&(t=e.length);for(var r=0,n=new Array(t);r<t;r++)n[r]=e[r];return n}},1106:function(e,t){e.exports=function(){throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}},1107:function(e,t,r){"use strict";(function(t){var n=r(997),a=n(r(998)),i=n(r(999)),u=n(r(1011)),s=n(r(1012)),o=r(680),c=r(1039),f=function(){function e(t){var r=t.filehandle,n=t.path;if((0,u.default)(this,e),r)this.filehandle=r;else{if(!n)throw new TypeError("either filehandle or path must be defined");this.filehandle=new c(n)}}return(0,s.default)(e,[{key:"_readLongWithOverflow",value:function(e){var t=arguments.length>1&&void 0!==arguments[1]?arguments[1]:0,r=!(arguments.length>2&&void 0!==arguments[2])||arguments[2],n=o.fromBytesLE(e.slice(t,t+8),r);if(n.greaterThan(Number.MAX_SAFE_INTEGER)||n.lessThan(Number.MIN_SAFE_INTEGER))throw new TypeError("integer overflow");return n.toNumber()}},{key:"_getIndex",value:function(){return this.index||(this.index=this._readIndex()),this.index}},{key:"_readIndex",value:function(){var e=(0,i.default)(a.default.mark((function e(){var r,n,i,u,s,o,c;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return r=t.allocUnsafe(8),e.next=3,this.filehandle.read(r,0,8,0);case 3:if(n=this._readLongWithOverflow(r,0,!0)){e.next=6;break}return e.abrupt("return",[[0,0]]);case 6:if((i=new Array(n+1))[0]=[0,0],!((u=16*n)>Number.MAX_SAFE_INTEGER)){e.next=11;break}throw new TypeError("integer overflow");case 11:return r=t.allocUnsafe(u),e.next=14,this.filehandle.read(r,0,u,8);case 14:for(s=0;s<n;s+=1)o=this._readLongWithOverflow(r,16*s),c=this._readLongWithOverflow(r,16*s+8),i[s+1]=[o,c];return e.abrupt("return",i);case 16:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"getLastBlock",value:function(){var e=(0,i.default)(a.default.mark((function e(){var t;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this._getIndex();case 2:if((t=e.sent).length){e.next=5;break}return e.abrupt("return",void 0);case 5:return e.abrupt("return",t[t.length-1]);case 6:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"getRelevantBlocksForRead",value:function(){var e=(0,i.default)(a.default.mark((function e(t,r){var n,i,u,s,o,c,f,l,h;return a.default.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:if(n=r+t,0!==t){e.next=3;break}return e.abrupt("return",[]);case 3:return e.next=5,this._getIndex();case 5:for(i=e.sent,u=[],s=function(e,t){var n=e[1],a=t?t[1]:1/0;return n<=r&&a>r?0:n<r?-1:1},o=0,c=i.length-1,f=Math.floor(i.length/2),l=s(i[f],i[f+1]);0!==l;)l>0?c=f-1:l<0&&(o=f+1),f=Math.ceil((c-o)/2)+o,l=s(i[f],i[f+1]);u.push(i[f]),h=f+1;case 15:if(!(h<i.length)){e.next=22;break}if(u.push(i[h]),!(i[h][1]>=n)){e.next=19;break}return e.abrupt("break",22);case 19:h+=1,e.next=15;break;case 22:return u[u.length-1][1]<n&&u.push([]),e.abrupt("return",u);case 24:case"end":return e.stop()}}),e,this)})));return function(t,r){return e.apply(this,arguments)}}()}]),e}();e.exports=f}).call(this,r(74).Buffer)},997:function(e,t){e.exports=function(e){return e&&e.__esModule?e:{default:e}}},998:function(e,t,r){e.exports=r(374)},999:function(e,t){function r(e,t,r,n,a,i,u){try{var s=e[i](u),o=s.value}catch(c){return void r(c)}s.done?t(o):Promise.resolve(o).then(n,a)}e.exports=function(e){return function(){var t=this,n=arguments;return new Promise((function(a,i){var u=e.apply(t,n);function s(e){r(u,a,i,s,o,"next",e)}function o(e){r(u,a,i,s,o,"throw",e)}s(void 0)}))}}}});
//# sourceMappingURL=15.3481a912c94ef75a9f78.worker.js.map