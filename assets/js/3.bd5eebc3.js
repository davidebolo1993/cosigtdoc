(window.webpackJsonp=window.webpackJsonp||[]).push([[3,21,24,31,33,34,35,44,47],{246:function(t,e,s){},248:function(t,e){t.exports=function(t){return null==t}},250:function(t,e,s){},251:function(t,e,s){},252:function(t,e,s){},253:function(t,e,s){},254:function(t,e,s){},261:function(t,e,s){"use strict";s(246)},262:function(t,e,s){"use strict";s.r(e);var i=s(288),a=s(268),n=s(242);function r(t,e){if("group"===e.type){const s=e.path&&Object(n.isActive)(t,e.path),i=e.children.some(e=>"group"===e.type?r(t,e):"page"===e.type&&Object(n.isActive)(t,e.path));return s||i}return!1}var o={name:"SidebarLinks",components:{SidebarGroup:i.default,SidebarLink:a.default},props:["items","depth","sidebarDepth","initialOpenGroupIndex"],data(){return{openGroupIndex:this.initialOpenGroupIndex||0}},watch:{$route(){this.refreshIndex()}},created(){this.refreshIndex()},methods:{refreshIndex(){const t=function(t,e){for(let s=0;s<e.length;s++){const i=e[s];if(r(t,i))return s}return-1}(this.$route,this.items);t>-1&&(this.openGroupIndex=t)},toggleGroup(t){this.openGroupIndex=t===this.openGroupIndex?-1:t},isActive(t){return Object(n.isActive)(this.$route,t.regularPath)}}},l=s(15),c=Object(l.a)(o,(function(){var t=this,e=t._self._c;return t.items.length?e("ul",{staticClass:"sidebar-links"},t._l(t.items,(function(s,i){return e("li",{key:i},["group"===s.type?e("SidebarGroup",{attrs:{item:s,open:i===t.openGroupIndex,collapsable:s.collapsable||s.collapsible,depth:t.depth},on:{toggle:function(e){return t.toggleGroup(i)}}}):e("SidebarLink",{attrs:{"sidebar-depth":t.sidebarDepth,item:s}})],1)})),0):t._e()}),[],!1,null,null,null);e.default=c.exports},263:function(t,e,s){"use strict";s.r(e);var i={name:"FluentWindowDevEdit20Regular"},a=s(15),n=Object(a.a)(i,(function(){var t=this._self._c;return t("svg",{attrs:{width:"1em",height:"1em",viewBox:"0 0 20 20"}},[t("g",{attrs:{fill:"none"}},[t("path",{attrs:{d:"M4.5 2A2.5 2.5 0 0 0 2 4.5v9A2.5 2.5 0 0 0 4.5 16h4.975c.11-.361.283-.7.51-1H4.5A1.5 1.5 0 0 1 3 13.5V5.999h12v3.944l.102-.102c.266-.267.572-.47.898-.61V4.5A2.5 2.5 0 0 0 13.5 2h-9zM2.999 4.499a1.5 1.5 0 0 1 1.5-1.5h9a1.5 1.5 0 0 1 1.5 1.5v.5h-12v-.5zm5.353 2.646a.5.5 0 0 1 0 .707L6.206 10l2.146 2.146a.5.5 0 0 1-.707.707l-2.5-2.5a.5.5 0 0 1 0-.707l2.5-2.5a.5.5 0 0 1 .707 0zm1.794 5a.5.5 0 0 0 .708.707l2.5-2.5a.5.5 0 0 0 0-.707l-2.5-2.5a.5.5 0 0 0-.708.707L12.293 10l-2.147 2.146zm5.663-1.597l-4.83 4.83a2.197 2.197 0 0 0-.577 1.02l-.375 1.498a.89.89 0 0 0 1.079 1.078l1.498-.374c.386-.097.739-.296 1.02-.578l4.83-4.83a1.87 1.87 0 0 0-2.645-2.644z",fill:"currentColor"}})])])}),[],!1,null,null,null);e.default=n.exports},267:function(t,e,s){"use strict";s.r(e);var i=s(15),a=Object(i.a)({},(function(){var t=this._self._c;return t("svg",{attrs:{width:"1em",height:"1em",viewBox:"0 0 24 24"}},[t("path",{attrs:{d:"M21 10.12h-6.78l2.74-2.82c-2.73-2.7-7.15-2.8-9.88-.1a6.887 6.887 0 0 0 0 9.8c2.73 2.7 7.15 2.7 9.88 0c1.36-1.35 2.04-2.92 2.04-4.9h2c0 1.98-.88 4.55-2.64 6.29c-3.51 3.48-9.21 3.48-12.72 0c-3.5-3.47-3.53-9.11-.02-12.58a8.987 8.987 0 0 1 12.65 0L21 3v7.12M12.5 8v4.25l3.5 2.08l-.72 1.21L11 13V8h1.5z",fill:"currentColor"}})])}),[],!1,null,null,null);e.default=a.exports},268:function(t,e,s){"use strict";s.r(e);var i=s(242);function a(t,e,s,i,a){const n={props:{to:e,activeClass:"",exactActiveClass:""},class:{active:i,"sidebar-link":!0}};return a>2&&(n.style={"padding-left":a+"rem"}),t("RouterLink",n,s)}function n(t,e,s,r,o,l=1){return!e||l>o?null:t("ul",{class:"sidebar-sub-headers"},e.map(e=>{const c=Object(i.isActive)(r,s+"#"+e.slug);return t("li",{class:"sidebar-sub-header"},[a(t,s+"#"+e.slug,e.title,c,e.level-1),n(t,e.children,s,r,o,l+1)])}))}var r={functional:!0,props:["item","sidebarDepth"],render(t,{parent:{$page:e,$site:s,$route:r,$themeConfig:o,$themeLocaleConfig:l},props:{item:c,sidebarDepth:u}}){const h=Object(i.isActive)(r,c.path),d="auto"===c.type?h||c.children.some(t=>Object(i.isActive)(r,c.basePath+"#"+t.slug)):h,p="external"===c.type?function(t,e,s){return t("a",{attrs:{href:e,target:"_blank",rel:"noopener noreferrer"},class:{"sidebar-link":!0}},[s,t("VPIconExternalLink")])}(t,c.path,c.title||c.path):a(t,c.path,c.title||c.path,d),f=[e.frontmatter.sidebarDepth,u,l.sidebarDepth,o.sidebarDepth,0].find(t=>void 0!==t),g=l.displayAllHeaders||o.displayAllHeaders;if("auto"===c.type)return[p,n(t,c.children,c.basePath,r,f)];if((d||g)&&c.headers&&!i.hashRE.test(c.path)){return[p,n(t,Object(i.groupHeaders)(c.headers),c.path,r,f)]}return p}},o=(s(261),s(15)),l=Object(o.a)(r,void 0,void 0,!1,null,null,null);e.default=l.exports},269:function(t,e,s){},272:function(t,e,s){},276:function(t,e,s){"use strict";s(250)},277:function(t,e,s){},278:function(t,e,s){"use strict";s(251)},279:function(t,e,s){var i=s(12),a=s(6),n=s(11);t.exports=function(t){return"string"==typeof t||!a(t)&&n(t)&&"[object String]"==i(t)}},280:function(t,e,s){"use strict";s(252)},281:function(t,e,s){},282:function(t,e,s){"use strict";s(253)},283:function(t,e,s){},284:function(t,e,s){},285:function(t,e,s){"use strict";s(254)},286:function(t,e,s){},288:function(t,e,s){"use strict";s.r(e);var i=s(242),a=s(255),n=s(256),r={name:"SidebarGroup",components:{DropdownTransition:a.default,Arrow:n.default},props:["item","open","collapsable","depth"],beforeCreate(){this.$options.components.SidebarLinks=s(262).default},methods:{isActive:i.isActive}},o=(s(282),s(15)),l=Object(o.a)(r,(function(){var t=this,e=t._self._c;return e("section",{staticClass:"sidebar-group",class:[{collapsable:t.collapsable,"is-sub-group":0!==t.depth},"depth-"+t.depth]},[t.item.path?e("RouterLink",{staticClass:"sidebar-heading clickable",class:{open:t.open,active:t.isActive(t.$route,t.item.path)},attrs:{to:t.item.path},nativeOn:{click:function(e){return t.$emit("toggle")}}},[e("span",[t._v(t._s(t.item.title))]),t._v(" "),t.collapsable?e("Arrow",{attrs:{open:t.open}}):t._e()],1):e("p",{staticClass:"sidebar-heading",class:{open:t.open},on:{click:function(e){return t.$emit("toggle")}}},[e("span",[t._v(t._s(t.item.title))]),t._v(" "),t.collapsable?e("Arrow",{attrs:{open:t.open}}):t._e()],1),t._v(" "),e("DropdownTransition",[t.open||!t.collapsable?e("SidebarLinks",{staticClass:"sidebar-group-items",attrs:{items:t.item.children,"sidebar-depth":t.item.sidebarDepth,"initial-open-group-index":t.item.initialOpenGroupIndex,depth:t.depth+1}}):t._e()],1)],1)}),[],!1,null,null,null);e.default=l.exports},291:function(t,e,s){"use strict";s.r(e);var i=s(248),a=s.n(i),n=s(242),r=s(263),o=s(267),l={name:"PageEdit",components:{VPIconEdit:r.default,VPIconLastUpdated:o.default},computed:{lastUpdated(){return this.$page.lastUpdated},lastUpdatedText(){return"string"==typeof this.$themeLocaleConfig.lastUpdated?this.$themeLocaleConfig.lastUpdated:"string"==typeof this.$site.themeConfig.lastUpdated?this.$site.themeConfig.lastUpdated:"Last Updated"},editLink(){const t=a()(this.$page.frontmatter.editLink)?this.$site.themeConfig.editLinks:this.$page.frontmatter.editLink,{repo:e,docsDir:s="",docsBranch:i="master",docsRepo:n=e}=this.$site.themeConfig;return t&&n&&this.$page.relativePath?this.createEditLink(e,n,s,i,this.$page.relativePath):null},editLinkText(){return this.$themeLocaleConfig.editLinkText||this.$site.themeConfig.editLinkText||"Edit this page"}},methods:{createEditLink(t,e,s,i,a){if(/bitbucket.org/.test(e)){return e.replace(n.endingSlashRE,"")+"/src"+`/${i}/`+(s?s.replace(n.endingSlashRE,"")+"/":"")+a+`?mode=edit&spa=0&at=${i}&fileviewer=file-view-default`}if(/gitlab.com/.test(e)){return e.replace(n.endingSlashRE,"")+"/-/edit"+`/${i}/`+(s?s.replace(n.endingSlashRE,"")+"/":"")+a}return(n.outboundRE.test(e)?e:"https://github.com/"+e).replace(n.endingSlashRE,"")+"/edit"+`/${i}/`+(s?s.replace(n.endingSlashRE,"")+"/":"")+a}}},c=(s(278),s(15)),u=Object(c.a)(l,(function(){var t=this,e=t._self._c;return e("footer",{staticClass:"page-edit"},[t.editLink?e("div",{staticClass:"edit-link"},[e("VPIconEdit",{staticClass:"edit-icon"}),t._v(" "),e("a",{attrs:{href:t.editLink,target:"_blank",rel:"noopener noreferrer"}},[t._v(t._s(t.editLinkText))])],1):t._e(),t._v(" "),t.lastUpdated?e("div",{staticClass:"last-updated"},[e("VPIconLastUpdated",{staticClass:"last-updated-icon"}),t._v(" "),e("span",{staticClass:"prefix"},[t._v(t._s(t.lastUpdatedText)+": ")]),t._v(" "),e("span",{staticClass:"time"},[t._v(t._s(t.lastUpdated))])],1):t._e()])}),[],!1,null,null,null);e.default=u.exports},292:function(t,e,s){"use strict";s.r(e);s(91);var i=s(242),a=s(279),n=s.n(a),r=s(248),o=s.n(r);const l={title:"N/A",path:"N/A",class:"noop"};var c={name:"PageNav",props:["sidebarItems"],computed:{prev(){return h(u.PREV,this)||l},next(){return h(u.NEXT,this)||l}}};const u={NEXT:{resolveLink:function(t,e){return d(t,e,1)},getThemeLinkConfig:({nextLinks:t})=>t,getPageLinkConfig:({frontmatter:t})=>t.next},PREV:{resolveLink:function(t,e){return d(t,e,-1)},getThemeLinkConfig:({prevLinks:t})=>t,getPageLinkConfig:({frontmatter:t})=>t.prev}};function h(t,{$themeConfig:e,$page:s,$route:a,$site:r,sidebarItems:l}){const{resolveLink:c,getThemeLinkConfig:u,getPageLinkConfig:h}=t,d=u(e),p=h(s),f=o()(p)?d:p;return!1===f?void 0:n()(f)?Object(i.resolvePage)(r.pages,f,a.path):c(s,l)}function d(t,e,s){const i=[];!function t(e,s){for(let i=0,a=e.length;i<a;i++)"group"===e[i].type?t(e[i].children||[],s):s.push(e[i])}(e,i);for(let e=0;e<i.length;e++){const a=i[e];if("page"===a.type&&a.path===decodeURIComponent(t.path))return i[e+s]}}var p=c,f=(s(280),s(15)),g=Object(f.a)(p,(function(){var t=this,e=t._self._c;return t.prev||t.next?e("div",{staticClass:"page-nav"},[e("div",{staticClass:"inner"},[t.prev?e("VPLink",{class:["prev",t.prev.class],attrs:{text:t.prev.title||t.prev.path,link:t.prev.path}},[e("span",{attrs:{slot:"before"},slot:"before"},[t._v("←")]),t._v(" "),e("span",{attrs:{slot:"after"},slot:"after"},[e("br"),t._v(" "),e("span",{staticClass:"prev-link"},[t._v(t._s(t.prev.path))])])]):t._e(),t._v(" "),t.next?e("VPLink",{class:["next",t.next.class],attrs:{text:t.next.title||t.next.path,link:t.next.path}},[e("span",{attrs:{slot:"after"},slot:"after"},[t._v("\n        →\n        "),e("br"),t._v(" "),e("span",{staticClass:"next-link"},[t._v(t._s(t.next.path))])])]):t._e()],1)]):t._e()}),[],!1,null,null,null);e.default=g.exports},293:function(t,e,s){"use strict";s.r(e);var i=s(242),a={props:["stick","tag"],data:()=>({needFloat:!1,stickBottom:0}),watch:{stick(){this.unStick(),this.stickHandle()}},methods:{stickHandle(){if(!this.stick)return;const t=Object(i.findContainerInVm)(this.stick,this);t&&(this._stickerScroll=()=>{const e=this.$el.getBoundingClientRect(),s=document.body.scrollTop+document.documentElement.scrollTop;this.needFloat=document.body.offsetHeight-s-e.height<t.offsetHeight,this.stickBottom=t.offsetHeight},this._stickerScroll(),window.addEventListener("scroll",this._stickerScroll))},unStick(){this.needFloat=!1,this.stickBottom=0,window.removeEventListener("scroll",this._stickerScroll)}},mounted(){this.stickHandle()},beforeDestroy(){this.unStick()}},n=(s(285),s(15)),r=Object(n.a)(a,(function(){return(0,this._self._c)(this.tag||"div",{tag:"component",staticClass:"sticker",class:this.needFloat?["stick-float"]:void 0,style:this.needFloat?{bottom:this.stickBottom+"px"}:void 0},[this._t("default")],2)}),[],!1,null,null,null);e.default=r.exports},294:function(t,e,s){"use strict";s.r(e);s(276);var i=s(15),a=Object(i.a)({},(function(){var t=this,e=t._self._c;return e("div",{staticClass:"sidebar-button",on:{click:function(e){return t.$emit("toggle-sidebar")}}},[e("svg",{staticClass:"icon",attrs:{xmlns:"http://www.w3.org/2000/svg","aria-hidden":"true",role:"img",viewBox:"0 0 448 512"}},[e("path",{attrs:{fill:"currentColor",d:"M436 124H12c-6.627 0-12-5.373-12-12V80c0-6.627 5.373-12 12-12h424c6.627 0 12 5.373 12 12v32c0 6.627-5.373 12-12 12zm0 160H12c-6.627 0-12-5.373-12-12v-32c0-6.627 5.373-12 12-12h424c6.627 0 12 5.373 12 12v32c0 6.627-5.373 12-12 12zm0 160H12c-6.627 0-12-5.373-12-12v-32c0-6.627 5.373-12 12-12h424c6.627 0 12 5.373 12 12v32c0 6.627-5.373 12-12 12z"}})])])}),[],!1,null,null,null);e.default=a.exports},295:function(t,e,s){"use strict";s(269)},301:function(t,e,s){"use strict";s(272)},303:function(t,e,s){"use strict";Object.defineProperty(e,"__esModule",{value:!0}),e.default={}},305:function(t,e,s){"use strict";s(277)},306:function(t,e,s){"use strict";s(281)},307:function(t,e,s){"use strict";s(283)},308:function(t,e,s){"use strict";s(284)},309:function(t,e,s){"use strict";s(286)},323:function(t,e,s){"use strict";s.r(e);var i={data:()=>({query:""}),created(){this.apiIndex=[]},computed:{apiGroups(){const t=this.$themeLocaleConfig.sidebar||this.$themeConfig.sidebar;if("object"==typeof t){const e=this.$localePath+"api/",s=t[e];if(s[0].path===e)return s.filter(t=>t.children).map(t=>({text:t.title,description:t.description,items:t.children.map(t=>{const e=this.$site.pages.find(e=>e.regularPath===t+".html"),{extractApiHeaders:s=[2,3]}=e.frontmatter;return{pageClass:e.frontmatter&&e.frontmatter.pageClass,text:e.title,link:e.path,headers:(e.headers||[]).filter(t=>s.includes(t.level))}})}))}return[]},filtered(){const t=this.query;return this.apiGroups.map(e=>{const s=e.items.map(({text:e,link:s,headers:i,pageClass:a})=>(i=i.filter(e=>e.slug.toLowerCase().includes(t.toLowerCase()))).length?{text:e,link:s,headers:i,pageClass:a}:null).filter(t=>t);return s.length?{text:e.text,id:e.text.replace(/\s+/g,"-"),description:e.description,pageClass:e.pageClass,items:s}:null}).filter(t=>t)}}},a=(s(295),s(15)),n=Object(a.a)(i,(function(){var t=this,e=t._self._c;return e("div",{attrs:{id:"api-index"}},[e("div",{staticClass:"header"},[e("h1",[t._v("API Reference")]),t._v(" "),e("input",{directives:[{name:"model",rawName:"v-model",value:t.query,expression:"query"}],staticClass:"api-filter",attrs:{placeholder:"Filter"},domProps:{value:t.query},on:{input:function(e){e.target.composing||(t.query=e.target.value)}}})]),t._v(" "),t._l(t.filtered,(function(s){return e("div",{key:s.id,staticClass:"api-section"},[e("h2",{attrs:{id:s.id}},[t._v("\n      "+t._s(s.text)+"\n      "),e("a",{staticClass:"header-anchor",attrs:{href:"#"+s.id,"aria-hidden":"true"}},[t._v("#")])]),t._v(" "),s.description?e("div",{staticClass:"tip custom-block api-section-description"},[e("p",{domProps:{innerHTML:t._s(s.description)}})]):t._e(),t._v(" "),e("div",{staticClass:"api-groups"},t._l(s.items,(function(s){return e("div",{key:s.text,staticClass:"api-group"},[e("h3",[t._v(t._s(s.text))]),t._v(" "),e("ul",{class:{"api-group-ul":!0,[s.pageClass]:!!s.pageClass}},t._l(s.headers,(function(t){return e("li",{key:t.slug,class:{"api-group-li":!0,["level-"+t.level]:!0}},[e("VPLink",{attrs:{text:t.title,link:s.link+"#"+t.slug}})],1)})),0)])})),0)])}))],2)}),[],!1,null,"31a87bd6",null);e.default=n.exports},325:function(t,e,s){"use strict";s.r(e);var i={computed:{data(){return this.$page.frontmatter}}},a=(s(301),s(15)),n=Object(a.a)(i,(function(){var t=this,e=t._self._c;return e("div",{staticClass:"home"},[e("section",{attrs:{id:"hero"}},[t.data.heroImage?e("img",{staticClass:"hero-img",attrs:{src:t.$withBase(t.data.heroImage),alt:t.data.heroAlt||"hero"}}):t._e(),t._v(" "),e("br"),t._v(" "),t.data.heroText?e("h1",{staticClass:"heroText"},[e("span",[t._v("\n        "+t._s(t.data.heroText||t.$description||"Welcome to your VuePress site")+"\n      ")])]):e("Content",{staticClass:"heroText",attrs:{"slot-key":"heroText"}}),t._v(" "),t.data.tagline?e("p",{staticClass:"tagline"},[t._v("\n      "+t._s(t.data.tagline)+"\n    ")]):e("Content",{staticClass:"tagline",attrs:{"slot-key":"tagline"}}),t._v(" "),e("p",{staticClass:"actions"},[t.data.actionText&&t.data.actionLink?e("VPLink",{staticClass:"action-link",attrs:{text:t.data.actionText,link:t.data.actionLink}},[e("svg",{staticClass:"icon",attrs:{slot:"after",xmlns:"http://www.w3.org/2000/svg",width:"10",height:"10",viewBox:"0 0 24 24"},slot:"after"},[e("path",{attrs:{d:"M13.025 1l-2.847 2.828 6.176 6.176h-16.354v3.992h16.354l-6.176 6.176 2.847 2.828 10.975-11z"}})])]):t._e(),t._v(" "),t.data.subActionText&&t.data.subActionLink?e("VPLink",{staticClass:"sub-action-link",attrs:{text:t.data.subActionText,link:t.data.subActionLink}}):t._e()],1)],1),t._v(" "),e("section",{directives:[{name:"show",rawName:"v-show",value:t.data.sponsors,expression:"data.sponsors"}],attrs:{id:"special-sponsor"}},[e("span",{staticClass:"special-sponsor-title"},[t._v(t._s(t.data.sponsorsText||"Special Sponsor"))]),t._v(" "),t._l(t.data.sponsors,(function(s){return e("span",{key:s.title,staticClass:"special-sponsor-item"},[e("span",[t._v(t._s(s.title))]),t._v(" "),e("a",{attrs:{href:s.link}},[e("img",{attrs:{src:s.img}})])])}))],2),t._v(" "),t.data.features?e("section",{staticClass:"vt-box-container",attrs:{id:"highlights"}},t._l(t.data.features,(function(s){return e("div",{key:s.title,staticClass:"vt-box"},[e("h3",[t._v(t._s(s.title))]),t._v(" "),e("p",[t._v(t._s(s.details))])])})),0):t._e(),t._v(" "),t.data.footer?e("div",{staticClass:"footer"},[t._v("\n    "+t._s(t.data.footer)+"\n  ")]):e("Content",{staticClass:"footer",attrs:{"slot-key":"footer"}})],1)}),[],!1,null,"4a10a1b0",null);e.default=n.exports},328:function(t,e,s){"use strict";s.r(e);var i=s(303),a=s.n(i),n=s(322),r=s(294),o=s(289);function l(t,e){return t.ownerDocument.defaultView.getComputedStyle(t,null)[e]}var c={name:"Navbar",components:{SidebarButton:r.default,NavLinks:o.default,SearchBox:n.default,AlgoliaSearchBox:a.a},props:["shouldShowStatusBar"],data:()=>({linksWrapMaxWidth:null}),computed:{algolia(){return this.$themeLocaleConfig.algolia||this.$site.themeConfig.algolia||{}},isAlgoliaSearch(){return this.algolia&&this.algolia.apiKey&&this.algolia.indexName}},mounted(){const t=parseInt(l(this.$el,"paddingLeft"))+parseInt(l(this.$el,"paddingRight")),e=()=>{document.documentElement.clientWidth<719?this.linksWrapMaxWidth=null:this.linksWrapMaxWidth=this.$el.offsetWidth-t-(this.$refs.siteName&&this.$refs.siteName.offsetWidth||0)};e(),window.addEventListener("resize",e,!1)}},u=(s(305),s(15)),h=Object(u.a)(c,(function(){var t=this,e=t._self._c;return e("header",{staticClass:"navbar",class:{"navbar-down":t.shouldShowStatusBar}},[e("div",{staticClass:"navbar-container"},[e("SidebarButton",{on:{"toggle-sidebar":function(e){return t.$emit("toggle-sidebar")}}}),t._v(" "),e("RouterLink",{staticClass:"home-link",attrs:{to:t.$localePath}},[t.$site.themeConfig.logo?e("img",{staticClass:"logo",attrs:{src:t.$withBase(t.$site.themeConfig.logo),alt:t.$siteTitle}}):t._e(),t._v(" "),t.$siteTitle?e("span",{ref:"siteName",staticClass:"site-name",class:{"can-hide":t.$site.themeConfig.logo}},[t._v(t._s(t.$siteTitle))]):t._e()]),t._v(" "),e("div",{staticClass:"links",style:t.linksWrapMaxWidth?{"max-width":t.linksWrapMaxWidth+"px"}:{}},[t.isAlgoliaSearch?e("AlgoliaSearchBox",{attrs:{options:t.algolia}}):!1!==t.$site.themeConfig.search&&!1!==t.$page.frontmatter.search?e("SearchBox"):t._e(),t._v(" "),e("NavLinks",{staticClass:"can-hide"})],1)],1)])}),[],!1,null,null,null);e.default=h.exports},330:function(t,e,s){"use strict";s.r(e);var i=s(291),a=s(292),n={components:{PageEdit:i.default,PageNav:a.default},props:["sidebarItems"],computed:{isHomepage(){return this.$route.path===this.$localeConfig.path},pageName(){return this.$route.path.replace(/\//g,"")}}},r=(s(306),s(15)),o=Object(r.a)(n,(function(){var t=this,e=t._self._c;return e("main",{class:["page",t.isHomepage?"homepage":"",t.pageName]},[t._t("top"),t._v(" "),e("Content",{staticClass:"theme-default-content vp-doc"}),t._v(" "),!1!==this.$page.frontmatter.pageEdit?e("PageEdit"):t._e(),t._v(" "),e("PageNav",t._b({},"PageNav",{sidebarItems:t.sidebarItems},!1)),t._v(" "),t._t("bottom")],2)}),[],!1,null,null,null);e.default=o.exports},331:function(t,e,s){"use strict";s.r(e);var i=s(262),a=s(289),n={name:"Sidebar",components:{SidebarLinks:i.default,NavLinks:a.default},props:["items"]},r=(s(307),s(15)),o=Object(r.a)(n,(function(){var t=this._self._c;return t("aside",{staticClass:"sidebar"},[t("NavLinks"),this._v(" "),this._t("top"),this._v(" "),t("SidebarLinks",{attrs:{depth:0,items:this.items}}),this._v(" "),this._t("bottom")],2)}),[],!1,null,null,null);e.default=o.exports},332:function(t,e,s){"use strict";s.r(e);var i={name:"StatusBar",computed:{statusText(){return this.$frontmatter&&this.$frontmatter.status||this.$themeLocaleConfig.status}}},a=(s(308),s(15)),n=Object(a.a)(i,(function(){var t=this._self._c;return t("div",{staticClass:"statusbar"},[this.statusText.startsWith("$")?t("Content",{attrs:{"slot-key":this.statusText.slice(1)}}):t("span",[this._v("\n    "+this._s(this.statusText)+"\n  ")])],1)}),[],!1,null,null,null);e.default=n.exports},333:function(t,e,s){"use strict";s.r(e);let i;function a(t){return t&&t.getBoundingClientRect?t.getBoundingClientRect().top+document.body.scrollTop+document.documentElement.scrollTop:0}var n={components:{Sticker:s(293).default},data:()=>({activeIndex:0}),computed:{visible(){return this.$frontmatter&&!1!==this.$frontmatter.toc&&!!(this.$page&&this.$page.headers&&this.$page.headers.length)}},watch:{activeIndex(){const t=(this.$refs.chairTocItem||[])[this.activeIndex];if(!t)return;const e=t.getBoundingClientRect(),s=this.$el.getBoundingClientRect(),i=e.top-s.top;i<20?this.$el.scrollTop=this.$el.scrollTop+i-20:i+e.height>s.height&&(this.$el.scrollTop+=e.top-(s.height-e.height))},$route(){}},methods:{onScroll(){void 0===i&&(i=a(this.$el));const t=document.body.scrollTop+document.documentElement.scrollTop,e=this.$page.headers||[];let s=0;const n=t=>{this.activeIndex=t};for(;s<e.length;s++){if(!(a(document.getElementById(e[s].slug))-50<t)){s||n(s);break}n(s)}},triggerEvt(){this._onScroll(),this._onHashChange()}},mounted(){const t=()=>{this.$emit("visible-change",this.visible)};t(),this.$watch("visible",t),setTimeout(()=>this.triggerEvt(),1e3),this._onScroll=()=>this.onScroll(),this._onHashChange=()=>{const t=decodeURIComponent(location.hash.substring(1)),e=(this.$page.headers||[]).findIndex(e=>e.slug===t);e>=0&&(this.activeIndex=e);const s=t&&document.getElementById(t);s&&window.scrollTo(0,a(s)-20)},window.addEventListener("scroll",this._onScroll)},beforeDestroy(){window.removeEventListener("scroll",this._onScroll),window.removeEventListener("hashchange",this._onHashChange)}},r=(s(309),s(15)),o=Object(r.a)(n,(function(){var t=this,e=t._self._c;return t.visible?e("Sticker",t._b({staticClass:"vuepress-toc"},"Sticker",t.$attrs,!1),[e("div",{staticClass:"on-this-page"},[t._v("ON THIS PAGE")]),t._v(" "),t._l(t.$page.headers,(function(s,i){return e("div",{key:i,ref:"chairTocItem",refInFor:!0,staticClass:"vuepress-toc-item",class:["vuepress-toc-h"+s.level,{active:t.activeIndex===i}]},[e("a",{attrs:{href:"#"+s.slug,title:s.title}},[t._v(t._s(s.title))])])}))],2):t._e()}),[],!1,null,null,null);e.default=o.exports},369:function(t,e,s){"use strict";s.r(e);var i=s(0),a=s(325),n=s(332),r=s(328),o=s(330),l=s(331),c=s(333),u=s(323),h=s(242),d={name:"Layout",components:{Home:a.default,Page:o.default,Sidebar:l.default,StatusBar:n.default,Navbar:r.default,Toc:c.default,API:u.default},data:()=>({isSidebarOpen:!1}),computed:{shouldShowStatusBar(){return!(!this.$frontmatter||!this.$frontmatter.status)||!!this.$themeLocaleConfig.status},shouldShowNavbar(){const{themeConfig:t}=this.$site,{frontmatter:e}=this.$page;return!1!==e.navbar&&!1!==t.navbar&&(this.$title||t.logo||t.repo||t.nav||this.$themeLocaleConfig.nav)},showShowSideSection(){const{frontmatter:t}=this.$page;return!t.home&&!t.api&&!t.pageLayout},shouldShowSidebar(){const{frontmatter:t}=this.$page;return this.showShowSideSection&&!1!==t.sidebar&&this.sidebarItems.length},shouldShowToc(){const{frontmatter:t}=this.$page;return this.showShowSideSection&&!1!==t.toc},sidebarItems(){return Object(h.resolveSidebarItems)(this.$page,this.$page.regularPath,this.$site,this.$localePath)},pageClasses(){const t=this.$page.frontmatter.pageClass;return[{"no-navbar":!this.shouldShowNavbar,"sidebar-open":this.isSidebarOpen,"no-sidebar":!this.shouldShowSidebar},t]},pageLayout(){const t=this.getPageLayout();if(t)return i.a.component(t)}},mounted(){this.$router.afterEach(()=>{this.isSidebarOpen=!1})},methods:{toggleSidebar(t){this.isSidebarOpen="boolean"==typeof t?t:!this.isSidebarOpen,this.$emit("toggle-sidebar",this.isSidebarOpen)},onTouchStart(t){this.touchStart={x:t.changedTouches[0].clientX,y:t.changedTouches[0].clientY}},onTouchEnd(t){const e=t.changedTouches[0].clientX-this.touchStart.x,s=t.changedTouches[0].clientY-this.touchStart.y;Math.abs(e)>Math.abs(s)&&Math.abs(e)>40&&(e>0&&this.touchStart.x<=80?this.toggleSidebar(!0):this.toggleSidebar(!1))},getPageLayout(){const t=this.$page.frontmatter.pageLayout;if(t&&(this.$vuepress.getLayoutAsyncComponent(t)||this.$vuepress.getVueComponent(t)))return t}}},p=s(15),f=Object(p.a)(d,(function(){var t=this,e=t._self._c;return e("div",{staticClass:"theme-container",class:t.pageClasses,on:{touchstart:t.onTouchStart,touchend:t.onTouchEnd}},[t.shouldShowStatusBar?e("StatusBar"):t._e(),t._v(" "),t.shouldShowNavbar?e("Navbar",{attrs:{shouldShowStatusBar:t.shouldShowStatusBar},on:{"toggle-sidebar":t.toggleSidebar}}):t._e(),t._v(" "),e("div",{staticClass:"sidebar-mask",on:{click:function(e){return t.toggleSidebar(!1)}}}),t._v(" "),e("Sidebar",{attrs:{items:t.sidebarItems},on:{"toggle-sidebar":t.toggleSidebar},scopedSlots:t._u([{key:"top",fn:function(){return[t._t("sidebar-top")]},proxy:!0},{key:"bottom",fn:function(){return[t._t("sidebar-bottom")]},proxy:!0}],null,!0)}),t._v(" "),t.$page.frontmatter.home?e("Home"):t.$page.frontmatter.api?e("API"):t.pageLayout?e("div",{staticClass:"custom-page"},[e(t.pageLayout,{tag:"component"})],1):e("Page",{attrs:{"sidebar-items":t.sidebarItems},scopedSlots:t._u([{key:"top",fn:function(){return[t._t("page-top")]},proxy:!0},{key:"bottom",fn:function(){return[t._t("page-bottom")]},proxy:!0}],null,!0)}),t._v(" "),t.shouldShowToc?e("Toc"):t._e()],1)}),[],!1,null,null,null);e.default=f.exports}}]);