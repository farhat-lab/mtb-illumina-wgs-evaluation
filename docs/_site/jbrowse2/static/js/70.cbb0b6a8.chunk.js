(this["webpackJsonp@jbrowse/web"]=this["webpackJsonp@jbrowse/web"]||[]).push([[70],{2515:function(e,t,n){"use strict";n.r(t),n.d(t,"default",(function(){return v}));var r=n(82),a=n(88),c=n(87),i=n.n(c),s=n(91),o=n(86),u=n(90),f=n(93),p=n(95),w=n(127),d=n(140),b=n(135),l=n(81),h=n(1078),v=function(e){Object(f.a)(n,e);var t=Object(p.a)(n);function n(){var e;Object(o.a)(this,n);for(var r=arguments.length,a=new Array(r),c=0;c<r;c++)a[c]=arguments[c];return(e=t.call.apply(t,[this].concat(a))).windowSize=1e3,e.windowDelta=1e3,e.gcMode="content",e}return Object(u.a)(n,[{key:"configure",value:function(){var e=Object(s.a)(i.a.mark((function e(){var t,n,r;return i.a.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return n=Object(l.readConfObject)(this.config,"sequenceAdapter"),e.next=3,null===(t=this.getSubAdapter)||void 0===t?void 0:t.call(this,n);case 3:if(r=e.sent){e.next=6;break}throw new Error("Error getting subadapter");case 6:return e.abrupt("return",r.dataAdapter);case 7:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"getRefNames",value:function(){var e=Object(s.a)(i.a.mark((function e(){var t;return i.a.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,this.configure();case 2:return t=e.sent,e.abrupt("return",t.getRefNames());case 4:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"getFeatures",value:function(e,t){var n=this;return this.windowSize=1e3,this.windowDelta=1e3,this.gcMode="content",Object(d.ObservableCreate)(function(){var c=Object(s.a)(i.a.mark((function c(s){var o,u,f,p,w,d,l,v,g,j,O,k,x,m,y,S,z,A;return i.a.wrap((function(c){for(;;)switch(c.prev=c.next){case 0:return c.next=2,n.configure();case 2:if(o=c.sent,u=1===n.windowSize?1:n.windowSize/2,f=1===n.windowSize,p=e.start,w=e.end,p=Math.max(0,p-u),!((w+=u)<0||p>w)){c.next=11;break}return s.complete(),c.abrupt("return");case 11:return d=o.getFeatures(Object(a.a)(Object(a.a)({},e),{},{start:p,end:w}),t),c.next=14,d.pipe(Object(h.a)()).toPromise();case 14:for(l=c.sent,v=Object(r.a)(l,1),g=v[0],j=g.get("seq"),O=u;O<j.length-u;O+=n.windowDelta){for(k=f?j[O]:j.slice(O-u,O+u),x=0,m=0,y=0,S=0;S<k.length;S++)"c"===k[S]||"C"===k[S]?x++:"g"!==k[S]&&"G"!==k[S]||m++,"N"!==k[S]&&y++;z=p,A=void 0,"content"===n.gcMode?A=(m+x)/(y||1):"skew"===n.gcMode&&(A=(m-x)/(m+x||1)),s.next(new b.a({uniqueId:"".concat(n.id,"_").concat(z+O),start:z+O,end:z+O+n.windowDelta,score:A}))}s.complete();case 20:case"end":return c.stop()}}),c)})));return function(e){return c.apply(this,arguments)}}())}},{key:"freeResources",value:function(){}}]),n}(w.BaseFeatureDataAdapter);v.capabilities=["hasLocalStats"]}}]);
//# sourceMappingURL=70.cbb0b6a8.chunk.js.map