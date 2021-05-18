this.webpackChunk([30],{1797:function(e,a,t){e.exports=function(){var e=[],a=[],t={},n={},s={};function r(e){return"string"===typeof e?new RegExp("^"+e+"$","i"):e}function o(e,a){return e===a?a:e===e.toLowerCase()?a.toLowerCase():e===e.toUpperCase()?a.toUpperCase():e[0]===e[0].toUpperCase()?a.charAt(0).toUpperCase()+a.substr(1).toLowerCase():a.toLowerCase()}function i(e,a){return e.replace(/\$(\d{1,2})/g,(function(e,t){return a[t]||""}))}function l(e,a){return e.replace(a[0],(function(t,n){var s=i(a[1],arguments);return o(""===t?e[n-1]:t,s)}))}function u(e,a,n){if(!e.length||t.hasOwnProperty(e))return a;for(var s=n.length;s--;){var r=n[s];if(r[0].test(a))return l(a,r)}return a}function c(e,a,t){return function(n){var s=n.toLowerCase();return a.hasOwnProperty(s)?o(n,s):e.hasOwnProperty(s)?o(n,e[s]):u(s,n,t)}}function m(e,a,t,n){return function(n){var s=n.toLowerCase();return!!a.hasOwnProperty(s)||!e.hasOwnProperty(s)&&u(s,s,t)===s}}function d(e,a,t){return(t?a+" ":"")+(1===a?d.singular(e):d.plural(e))}return d.plural=c(s,n,e),d.isPlural=m(s,n,e),d.singular=c(n,s,a),d.isSingular=m(n,s,a),d.addPluralRule=function(a,t){e.push([r(a),t])},d.addSingularRule=function(e,t){a.push([r(e),t])},d.addUncountableRule=function(e){"string"!==typeof e?(d.addPluralRule(e,"$0"),d.addSingularRule(e,"$0")):t[e.toLowerCase()]=!0},d.addIrregularRule=function(e,a){a=a.toLowerCase(),e=e.toLowerCase(),s[e]=a,n[a]=e},[["I","we"],["me","us"],["he","they"],["she","they"],["them","them"],["myself","ourselves"],["yourself","yourselves"],["itself","themselves"],["herself","themselves"],["himself","themselves"],["themself","themselves"],["is","are"],["was","were"],["has","have"],["this","these"],["that","those"],["echo","echoes"],["dingo","dingoes"],["volcano","volcanoes"],["tornado","tornadoes"],["torpedo","torpedoes"],["genus","genera"],["viscus","viscera"],["stigma","stigmata"],["stoma","stomata"],["dogma","dogmata"],["lemma","lemmata"],["schema","schemata"],["anathema","anathemata"],["ox","oxen"],["axe","axes"],["die","dice"],["yes","yeses"],["foot","feet"],["eave","eaves"],["goose","geese"],["tooth","teeth"],["quiz","quizzes"],["human","humans"],["proof","proofs"],["carve","carves"],["valve","valves"],["looey","looies"],["thief","thieves"],["groove","grooves"],["pickaxe","pickaxes"],["passerby","passersby"]].forEach((function(e){return d.addIrregularRule(e[0],e[1])})),[[/s?$/i,"s"],[/[^\u0000-\u007F]$/i,"$0"],[/([^aeiou]ese)$/i,"$1"],[/(ax|test)is$/i,"$1es"],[/(alias|[^aou]us|t[lm]as|gas|ris)$/i,"$1es"],[/(e[mn]u)s?$/i,"$1s"],[/([^l]ias|[aeiou]las|[ejzr]as|[iu]am)$/i,"$1"],[/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i,"$1i"],[/(alumn|alg|vertebr)(?:a|ae)$/i,"$1ae"],[/(seraph|cherub)(?:im)?$/i,"$1im"],[/(her|at|gr)o$/i,"$1oes"],[/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|automat|quor)(?:a|um)$/i,"$1a"],[/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)(?:a|on)$/i,"$1a"],[/sis$/i,"ses"],[/(?:(kni|wi|li)fe|(ar|l|ea|eo|oa|hoo)f)$/i,"$1$2ves"],[/([^aeiouy]|qu)y$/i,"$1ies"],[/([^ch][ieo][ln])ey$/i,"$1ies"],[/(x|ch|ss|sh|zz)$/i,"$1es"],[/(matr|cod|mur|sil|vert|ind|append)(?:ix|ex)$/i,"$1ices"],[/\b((?:tit)?m|l)(?:ice|ouse)$/i,"$1ice"],[/(pe)(?:rson|ople)$/i,"$1ople"],[/(child)(?:ren)?$/i,"$1ren"],[/eaux$/i,"$0"],[/m[ae]n$/i,"men"],["thou","you"]].forEach((function(e){return d.addPluralRule(e[0],e[1])})),[[/s$/i,""],[/(ss)$/i,"$1"],[/(wi|kni|(?:after|half|high|low|mid|non|night|[^\w]|^)li)ves$/i,"$1fe"],[/(ar|(?:wo|[ae])l|[eo][ao])ves$/i,"$1f"],[/ies$/i,"y"],[/\b([pl]|zomb|(?:neck|cross)?t|coll|faer|food|gen|goon|group|lass|talk|goal|cut)ies$/i,"$1ie"],[/\b(mon|smil)ies$/i,"$1ey"],[/\b((?:tit)?m|l)ice$/i,"$1ouse"],[/(seraph|cherub)im$/i,"$1"],[/(x|ch|ss|sh|zz|tto|go|cho|alias|[^aou]us|t[lm]as|gas|(?:her|at|gr)o|[aeiou]ris)(?:es)?$/i,"$1"],[/(analy|diagno|parenthe|progno|synop|the|empha|cri|ne)(?:sis|ses)$/i,"$1sis"],[/(movie|twelve|abuse|e[mn]u)s$/i,"$1"],[/(test)(?:is|es)$/i,"$1is"],[/(alumn|syllab|vir|radi|nucle|fung|cact|stimul|termin|bacill|foc|uter|loc|strat)(?:us|i)$/i,"$1us"],[/(agend|addend|millenni|dat|extrem|bacteri|desiderat|strat|candelabr|errat|ov|symposi|curricul|quor)a$/i,"$1um"],[/(apheli|hyperbat|periheli|asyndet|noumen|phenomen|criteri|organ|prolegomen|hedr|automat)a$/i,"$1on"],[/(alumn|alg|vertebr)ae$/i,"$1a"],[/(cod|mur|sil|vert|ind)ices$/i,"$1ex"],[/(matr|append)ices$/i,"$1ix"],[/(pe)(rson|ople)$/i,"$1rson"],[/(child)ren$/i,"$1"],[/(eau)x?$/i,"$1"],[/men$/i,"man"]].forEach((function(e){return d.addSingularRule(e[0],e[1])})),["adulthood","advice","agenda","aid","aircraft","alcohol","ammo","analytics","anime","athletics","audio","bison","blood","bream","buffalo","butter","carp","cash","chassis","chess","clothing","cod","commerce","cooperation","corps","debris","diabetes","digestion","elk","energy","equipment","excretion","expertise","firmware","flounder","fun","gallows","garbage","graffiti","hardware","headquarters","health","herpes","highjinks","homework","housework","information","jeans","justice","kudos","labour","literature","machinery","mackerel","mail","media","mews","moose","music","mud","manga","news","only","personnel","pike","plankton","pliers","police","pollution","premises","rain","research","rice","salmon","scissors","series","sewage","shambles","shrimp","software","species","staff","swine","tennis","traffic","transportation","trout","tuna","wealth","welfare","whiting","wildebeest","wildlife","you",/pok[e\xe9]mon$/i,/[^aeiou]ese$/i,/deer$/i,/fish$/i,/measles$/i,/o[iu]s$/i,/pox$/i,/sheep$/i].forEach(d.addUncountableRule),d}()},2154:function(e,a,t){"use strict";t.r(a),t.d(a,"default",(function(){return R}));var n=t(9),s=t(0),r=t.n(s),o=t(18),i=t(936),l=t(108),u=t(937),c=t(938),m=t(941),d=t(171),h=t(896),f=t(897),$=t(674),p=t(904),g=t(903),v=t(149),b=t(95),y=t(965),w=t(1797),E=t.n(w),x=t(458),C=Object(l.a)((function(e){return{root:{margin:e.spacing(1)},message:{padding:e.spacing(3)},titleBox:{color:"#fff",backgroundColor:e.palette.primary.main,textAlign:"center"},dialogContent:{width:600},resetButton:{justifyContent:"center",marginBottom:"6px"}}})),k=Object(o.observer)((function(e){var a=e.session,t=e.selectedDefault,n=e.handleRadio,s=C();return r.a.createElement(v.a,{className:s.root},r.a.createElement(h.a,{subheader:r.a.createElement(g.a,null,"Currently open session")},r.a.createElement(f.a,null,r.a.createElement($.a,null,r.a.createElement(y.a,{checked:a.name===t,onChange:function(){return n(a)}})),r.a.createElement(p.a,{primary:a.name}))))})),R=Object(o.observer)((function(e){var a=e.rootModel,t=e.open,o=e.onClose,l=e.currentDefault,w=C(),R=a.session,S=Object(s.useState)(l),j=Object(n.a)(S,2),O=j[0],z=j[1];function L(e){z(e.name),a.jbrowse.setDefaultSessionConf(e),R.notify("Set default session to ".concat(e.name),"success")}return r.a.createElement(i.a,{open:t},r.a.createElement(u.a,{className:w.titleBox},"Set Default Session"),r.a.createElement(c.a,null,r.a.createElement(x.a,{className:w.resetButton,container:!0},r.a.createElement(x.a,{item:!0},r.a.createElement(d.a,{color:"secondary",variant:"contained",onClick:function(){z("New session"),a.jbrowse.setDefaultSessionConf({name:"New session"}),R.notify("Reset default session","success")}},"Clear default session"))),r.a.createElement(k,{session:R,selectedDefault:O,handleRadio:L}),r.a.createElement(v.a,{className:w.root},r.a.createElement(h.a,{subheader:r.a.createElement(g.a,null,"Saved sessions")},R.savedSessions.length?R.savedSessions.map((function(e){var a=e.views,t=void 0===a?[]:a,n=t.map((function(e){return e.tracks.length})).reduce((function(e,a){return e+a}),0);return e.name!==R.name?r.a.createElement(f.a,{key:e.name},r.a.createElement($.a,null,r.a.createElement(y.a,{checked:e.name===O,onChange:function(){return L(e)}})),r.a.createElement(p.a,{primary:e.name,secondary:"".concat(t.length," ").concat(E()("view",t.length),"; ").concat(n,"\n                             open ").concat(E()("track",n))})):null})):r.a.createElement(b.a,{className:w.message},"No saved sessions found")))),r.a.createElement(m.a,null,r.a.createElement(d.a,{color:"secondary",variant:"contained",onClick:function(){o(!1)}},"Return")))}))}});
//# sourceMappingURL=30.3481a912c94ef75a9f78.worker.js.map