(this["webpackJsonp@jbrowse/web"]=this["webpackJsonp@jbrowse/web"]||[]).push([[38],{2543:function(e,a,t){"use strict";t.r(a),t.d(a,"default",(function(){return R}));var n=t(82),s=t(1),r=t.n(s),o=t(85),i=t(2454),l=t(51),u=t(2457),c=t(2458),m=t(2463),h=t(2437),d=t(2425),f=t(2432),$=t(2438),p=t(2451),g=t(2450),v=t(2423),b=t(2435),w=t(2532),y=t(850),E=t.n(y),x=t(2460),k=Object(l.a)((function(e){return{root:{margin:e.spacing(1)},message:{padding:e.spacing(3)},titleBox:{color:"#fff",backgroundColor:e.palette.primary.main,textAlign:"center"},dialogContent:{width:600},resetButton:{justifyContent:"center",marginBottom:"6px"}}})),C=Object(o.observer)((function(e){var a=e.session,t=e.selectedDefault,n=e.handleRadio,s=k();return r.a.createElement(v.a,{className:s.root},r.a.createElement(d.a,{subheader:r.a.createElement(g.a,null,"Currently open session")},r.a.createElement(f.a,null,r.a.createElement($.a,null,r.a.createElement(w.a,{checked:a.name===t,onChange:function(){return n(a)}})),r.a.createElement(p.a,{primary:a.name}))))})),R=Object(o.observer)((function(e){var a=e.rootModel,t=e.open,o=e.onClose,l=e.currentDefault,y=k(),R=a.session,j=Object(s.useState)(l),S=Object(n.a)(j,2),O=S[0],z=S[1];function L(e){z(e.name),a.jbrowse.setDefaultSessionConf(e),R.notify("Set default session to ".concat(e.name),"success")}return r.a.createElement(i.a,{open:t},r.a.createElement(u.a,{className:y.titleBox},"Set Default Session"),r.a.createElement(c.a,null,r.a.createElement(x.a,{className:y.resetButton,container:!0},r.a.createElement(x.a,{item:!0},r.a.createElement(h.a,{color:"secondary",variant:"contained",onClick:function(){z("New session"),a.jbrowse.setDefaultSessionConf({name:"New session"}),R.notify("Reset default session","success")}},"Clear default session"))),r.a.createElement(C,{session:R,selectedDefault:O,handleRadio:L}),r.a.createElement(v.a,{className:y.root},r.a.createElement(d.a,{subheader:r.a.createElement(g.a,null,"Saved sessions")},R.savedSessions.length?R.savedSessions.map((function(e){var a=e.views,t=void 0===a?[]:a,n=t.map((function(e){return e.tracks.length})).reduce((function(e,a){return e+a}),0);return e.name!==R.name?r.a.createElement(f.a,{key:e.name},r.a.createElement($.a,null,r.a.createElement(w.a,{checked:e.name===O,onChange:function(){return L(e)}})),r.a.createElement(p.a,{primary:e.name,secondary:"".concat(t.length," ").concat(E()("view",t.length),"; ").concat(n,"\n                             open ").concat(E()("track",n))})):null})):r.a.createElement(b.a,{className:y.message},"No saved sessions found")))),r.a.createElement(m.a,null,r.a.createElement(h.a,{color:"secondary",variant:"contained",onClick:function(){o(!1)}},"Return")))}))},850:function(e,a,t){e.exports=function(){var e=[],a=[],t={},n={},s={};function r(e){return"string"===typeof e?new RegExp("^"+e+"$","i"):e}function o(e,a){return e===a?a:e===e.toLowerCase()?a.toLowerCase():e===e.toUpperCase()?a.toUpperCase():e[0]===e[0].toUpperCase()?a.charAt(0).toUpperCase()+a.substr(1).toLowerCase():a.toLowerCase()}function i(e,a){return e.replace(/\$(\d{1,2})/g,(function(e,t){return a[t]||""}))}function l(e,a){return e.replace(a[0],(function(t,n){var s=i(a[1],arguments);return o(""===t?e[n-1]:t,s)}))}function u(e,a,n){if(!e.length||t.hasOwnProperty(e))return a;for(var s=n.length;s--;){var r=n[s];if(r[0].test(a))return l(a,r)}return a}function c(e,a,t){return function(n){var s=n.toLowerCase();return a.hasOwnProperty(s)?o(n,s):e.hasOwnProperty(s)?o(n,e[s]):u(s,n,t)}}function m(e,a,t,n){return function(n){var s=n.toLowerCase();return!!a.hasOwnProperty(s)||!e.hasOwnProperty(s)&&u(s,s,t)===s}}function h(e,a,t){return(t?a+" ":"")+(1===a?h.singular(e):h.plural(e))}return h.plural=c(s,n,e),h.isPlural=m(s,n,e),h.singular=c(n,s,a),h.isSingular=m(n,s,a),h.addPluralRule=function(a,t){e.push([r(a),t])},h.addSingularRule=function(e,t){a.push([r(e),t])},h.addUncountableRule=function(e){"string"!==typeof e?(h.addPluralRule(e,"$0"),h.addSingularRule(e,"$0")):t[e.toLowerCase()]=!0},h.addIrregularRule=function(e,a){a=a.toLowerCase(),e=e.toLowerCase(),s[e]=a,n[a]=e},[["I","we"],["me","us"],["he","they"],["she","they"],["them","them"],["myself","ourselves"],["yourself","yourselves"],["itself","themselves"],["herself","themselves"],["himself","themselves"],["themself","themselves"],["is","are"],["was","were"],["has","have"],["this","these"],["that","those"],["echo","echoes"],["dingo","dingoes"],["volcano","volcanoes"],["tornado","tornadoes"],["torpedo","torpedoes"],["genus","genera"],["viscus","viscera"],["stigma","stigmata"],["stoma","stomata"],["dogma","dogmata"],["lemma","lemmata"],["schema","schemata"],["anathema","anathemata"],["ox","oxen"],["axe","axes"],["die","dice"],["yes","yeses"],["foot","feet"],["eave","eaves"],["goose","geese"],["tooth","teeth"],["quiz","quizzes"],["human","humans"],["proof","proofs"],["carve","carves"],["valve","valves"],["looey","looies"],["thief","thieves"],["groove","grooves"],["pickaxe","pickaxes"],["passerby","passersby"]].forEach((function(e){return h.addIrregularRule(e[0],e[1])})),[[/s?$/i,"s"],[/[^\u0000-\u007F]$/i,"$0"],[/([^aeiou]ese)$/i,"$1"],[/(ax|test)is$/i,"$1es"],[/(alias|[^aou]us|t[lm]as|gas|ris)$/i,"$1es"],[/(e[mn]u)s?$/i,"$1s"],[/([^l]ias|[aeiou]las|[ejzr]as|[iu]am)$/i,"$1"],[/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i,"$1i"],[/(alumn|alg|vertebr)(?:a|ae)$/i,"$1ae"],[/(seraph|cherub)(?:im)?$/i,"$1im"],[/(her|at|gr)o$/i,"$1oes"],[/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|automat|quor)(?:a|um)$/i,"$1a"],[/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)(?:a|on)$/i,"$1a"],[/sis$/i,"ses"],[/(?:(kni|wi|li)fe|(ar|l|ea|eo|oa|hoo)f)$/i,"$1$2ves"],[/([^aeiouy]|qu)y$/i,"$1ies"],[/([^ch][ieo][ln])ey$/i,"$1ies"],[/(x|ch|ss|sh|zz)$/i,"$1es"],[/(matr|cod|mur|sil|vert|ind|append)(?:ix|ex)$/i,"$1ices"],[/\b((?:tit)?m|l)(?:ice|ouse)$/i,"$1ice"],[/(pe)(?:rson|ople)$/i,"$1ople"],[/(child)(?:ren)?$/i,"$1ren"],[/eaux$/i,"$0"],[/m[ae]n$/i,"men"],["thou","you"]].forEach((function(e){return h.addPluralRule(e[0],e[1])})),[[/s$/i,""],[/(ss)$/i,"$1"],[/(wi|kni|(?:after|half|high|low|mid|non|night|[^\w]|^)li)ves$/i,"$1fe"],[/(ar|(?:wo|[ae])l|[eo][ao])ves$/i,"$1f"],[/ies$/i,"y"],[/\b([pl]|zomb|(?:neck|cross)?t|coll|faer|food|gen|goon|group|lass|talk|goal|cut)ies$/i,"$1ie"],[/\b(mon|smil)ies$/i,"$1ey"],[/\b((?:tit)?m|l)ice$/i,"$1ouse"],[/(seraph|cherub)im$/i,"$1"],[/(x|ch|ss|sh|zz|tto|go|cho|alias|[^aou]us|t[lm]as|gas|(?:her|at|gr)o|[aeiou]ris)(?:es)?$/i,"$1"],[/(analy|diagno|parenthe|progno|synop|the|empha|cri|ne)(?:sis|ses)$/i,"$1sis"],[/(movie|twelve|abuse|e[mn]u)s$/i,"$1"],[/(test)(?:is|es)$/i,"$1is"],[/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i,"$1us"],[/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|quor)a$/i,"$1um"],[/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)a$/i,"$1on"],[/(alumn|alg|vertebr)ae$/i,"$1a"],[/(cod|mur|sil|vert|ind)ices$/i,"$1ex"],[/(matr|append)ices$/i,"$1ix"],[/(pe)(rson|ople)$/i,"$1rson"],[/(child)ren$/i,"$1"],[/(eau)x?$/i,"$1"],[/men$/i,"man"]].forEach((function(e){return h.addSingularRule(e[0],e[1])})),["adulthood","advice","agenda","aid","aircraft","alcohol","ammo","analytics","anime","athletics","audio","bison","blood","bream","buffalo","butter","carp","cash","chassis","chess","clothing","cod","commerce","cooperation","corps","debris","diabetes","digestion","elk","energy","equipment","excretion","expertise","firmware","flounder","fun","gallows","garbage","graffiti","hardware","headquarters","health","herpes","highjinks","homework","housework","information","jeans","justice","kudos","labour","literature","machinery","mackerel","mail","media","mews","moose","music","mud","manga","news","only","personnel","pike","plankton","pliers","police","pollution","premises","rain","research","rice","salmon","scissors","series","sewage","shambles","shrimp","software","species","staff","swine","tennis","traffic","transportation","trout","tuna","wealth","welfare","whiting","wildebeest","wildlife","you",/pok[e\xe9]mon$/i,/[^aeiou]ese$/i,/deer$/i,/fish$/i,/measles$/i,/o[iu]s$/i,/pox$/i,/sheep$/i].forEach(h.addUncountableRule),h}()}}]);
//# sourceMappingURL=38.731e9781.chunk.js.map