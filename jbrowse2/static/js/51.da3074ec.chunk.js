(this["webpackJsonp@jbrowse/web"]=this["webpackJsonp@jbrowse/web"]||[]).push([[51],{2498:function(e,t,a){"use strict";a.r(t);var o=a(82),r=a(1),n=a.n(r),l=a(85),c=a(51),s=a(2454),i=a(2457),m=a(2427),p=a(2458),u=a(2435),d=a(2412),g=a(2437),b=a(150),h=a.n(b),f=Object(c.a)((function(e){return{root:{width:300},closeButton:{position:"absolute",right:e.spacing(1),top:e.spacing(1),color:e.palette.grey[500]}}}));t.default=Object(l.observer)((function(e){var t=f(),a=e.model,l=e.handleClose,c=Object(r.useState)(""),b=Object(o.a)(c,2),E=b[0],v=b[1],w=E.match(/^[A-Za-z][A-Za-z0-9]$/);return n.a.createElement(s.a,{open:!0,onClose:l},n.a.createElement(i.a,null,"Color by tag",n.a.createElement(m.a,{"aria-label":"close",className:t.closeButton,onClick:l},n.a.createElement(h.a,null))),n.a.createElement(p.a,{style:{overflowX:"hidden"}},n.a.createElement("div",{className:t.root},n.a.createElement(u.a,null,"Enter tag to color by: "),n.a.createElement(u.a,{color:"textSecondary"},"Examples: XS or TS for RNA-seq inferred read strand, ts (lower-case) for minimap2 read strand, HP for haplotype, RG for read group, etc."),n.a.createElement(d.a,{value:E,onChange:function(e){v(e.target.value)},placeholder:"Enter tag name",inputProps:{maxLength:2,"data-testid":"color-tag-name-input"},error:2===E.length&&!w,helperText:2!==E.length||w?"":"Not a valid tag",autoComplete:"off","data-testid":"color-tag-name"}),n.a.createElement(g.a,{variant:"contained",color:"primary",style:{marginLeft:20},onClick:function(){a.setColorScheme({type:"tag",tag:E}),l()},disabled:!w},"Submit"))))}))}}]);
//# sourceMappingURL=51.da3074ec.chunk.js.map