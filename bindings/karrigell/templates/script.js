function get(query) {
   z = ($.ajax({
     url: 'api.ks/'+query,
     dataType: 'json',
     async: false
   }));

   return $.parseJSON(z.responseText);
}

function map(elements, f) {
    var res = [];
    for (var i=0; i<elements.length; i++)
        res.push(f(elements[i], i));
    return res;
}

function CachedMap(getter) {
   var self = this;
   self._storage = {};
   self.at = function (key) {
        if (self._storage[key] == undefined){
            //console.log(key + ' not found -> fetching');
            self._storage[key] = ko.observable(getter(key));
        } //else console.log(key + 'found');
        return self._storage[key];
   };  
   self.set = function (key, value) {
        self._storage[key](value);
   }
}

function KsCachedMap(query_fun) {
    return new CachedMap(function (key) {
        return get(query_fun(key));
    });
}

function ParamCachedMap(parent, entity, query_fun) {
    return new CachedMap(function (key) {
            return loadParams(get(query_fun(key)), parent, [entity]);
        });
}

function loadResult(results) {
    return map(results, function (r) {
        return { 
            label: ko.observable(r[0]),
            kind: r[1], 
            visible: ko.observable(r[2])
        }; 
    });
}

function ResultCachedMap(query_fun) {
    return new CachedMap(function (key) {
            return loadResult(get(query_fun(key)));
        });
}

var iterables = ko.observableArray([]);

if (!String.prototype.trim) {
    String.prototype.trim=function(){return this.replace(/^\s\s*/, '').replace(/\s\s*$/, '');};
}

function isInteger (s) {
    var isInteger_re     = /^\s*(\+|-)?\d+\s*$/;
    return String(s).search (isInteger_re) != -1
}

function isFloat (s) {
    var isDecimal_re     = /^\s*(\+|-)?((\d+(\.\d+)?)|(\.\d+))\s*$/;
    return String(s).search (isDecimal_re) != -1
}

function _parseInt(x) {
    if (typeof(x) == "string" && !isInteger(x.trim()))
        return NaN;
    return parseInt(x,10);
}

function _parseFloat(x) {
    if (typeof(x) == "string" && !isFloat(x.trim()))
        return NaN;
    return parseFloat(x);
}

function less(y, strict) {
    return strict ? (function (x) {
        return x < y ? x : NaN;
    }) : (function (x) {
        return x <= y ? x : NaN;
    });
}

function greater(y, strict) {
    return strict ? (function (x) {
        return x > y ? x : NaN;
    }) : (function (x) {
        return x >= y ? x : NaN;
    });
} 

function combine(f,g) {
    return function (x) {
        return f(g(x));
    }
}

function NaN2error(f) {
    return function (x) {
        var r = f(x);
        return !isNaN(r) ? r : "_error_";
    }
}

var lastLabel = ko.observable(undefined);
var lastValue = ko.observable(undefined);

function ScalarValue(value, label, converter, iterable, reloader) {
    var converter = converter || _parseFloat;
    var iterable = iterable || true;
    var reloader = reloader || null;
    var self = this;
    var nanconverter = NaN2error(converter);
    self.valueRaw = ko.observable(value);
    self.value = ko.computed({
        read: function () { return self.valueRaw(); },
        write: function (x) { 
            if (reloader) 
                reloader(nanconverter(x));
            else {
                self.valueRaw(x); 
                lastLabel(label);
                lastValue(x);
            }
        }
    })
    self.valueInvalid = ko.computed(function(){
        return isNaN(converter(self.value()));
    })

    self.iterable = iterable;

    self.iterateTo = ko.observable(value);
    self.iterationsCount = ko.observable(10);

    self.iterateToInvalid = ko.computed(function(){
        return isNaN(converter(self.iterateTo()));
    })

    var iterateToChecker = combine(greater(1, true),_parseInt);

    self.iterationCountInvalid = ko.computed(function(){
        return isNaN(iterateToChecker(self.iterationsCount()));
    })


    self.hasIteration = ko.computed({
        read: function () { return iterables().indexOf(self) != -1; },
        write: function (value) {
            if (value) {
                iterables.push(self);
                if (iterables().length > 2)
                    iterables.splice(0,1);
            } else {
                iterables.splice(iterables.indexOf(self), 1);
            }
        }
    })

    self.serialized = ko.computed(function() {
        var x = nanconverter(self.value());
        var to = nanconverter(self.iterateTo()); 
        return (self.hasIteration() && x != to
            ? ['iterate', 
                x, 
                to, 
                NaN2error(iterateToChecker)(self.iterationsCount())
            ] : x);
    });

    self.asHistory = ko.computed(function () {
        return (self.hasIteration() 
            ? "{" + self.value() + ":" + self.iterateTo() + "/" + self.iterationsCount() + "}" 
            : self.value());
    });
}

function renderConstraint(constraint) {
    var R = '&#8477;';
    var Z = '&#8484;';
    var and = '&#8745;';
    var inf = '&infin;';

    var type = constraint[0] == 'Z' ? Z : R;

    if (constraint.length == 3) {
        if (constraint[1].length==0 && constraint[2].length==0)
            return type;
        var lo = (constraint[1].length == 2 ? (constraint[1][0] == 'I' ? '[' : '(')+constraint[1][1] : ("(-"+inf));
        var hi = (constraint[2].length == 2 ? constraint[2][1]+(constraint[2][0] == 'I' ? ']' : ')') : ("+"+inf+")"));
        return type + and + lo + "," + hi;
    }
    return type;
}

function makeConverter(constraint) {

    var conv = constraint[0] == 'Z' ? _parseInt : _parseFloat;

    if (constraint.length == 3) {
        if (constraint[1].length == 2) 
            conv = combine(greater(constraint[1][1], constraint[1][0] == 'S'), conv);
        if (constraint[2].length == 2) 
            conv = combine(less(constraint[2][1], constraint[2][0] == 'S'), conv);
    }
    return conv;
}

function ScalarField(label, value, constraint, iterable, root, prefix, setter) {
    var root = root || null;
    var prefix = prefix || [];
    var setter = setter || "";
    var self = this;
    self.type = 0;
    self.label = label;
    var reloader = setter=="" ? null : function (x) {
        var c = prefix.slice(0);
        c.push(setter);
        root.reload(c, x);
    };
    self.hasTrivialSetter = setter.length == 0;
    self.value = new ScalarValue(value, label, makeConverter(constraint), iterable, reloader);
    self.renderer = 'scalar-row-template';

    self.getFields = function() {
        return [self];
    }

    self.serialized = function() {
        return [self.value.serialized()];
    }

    self.asHistory = function() {
        return [[self.label, self.value.asHistory()]];
    }

    self.constraint = renderConstraint(constraint);

}

function VectorField(label, values) {
    var self = this;
    self.label = label;
    self.type = 1;
    self.renderer = 'vector-row-template';

    self.values = $.map(values, function (value, i) { return new ScalarValue(value, label+'['+i+']'); });

    self.spansize = ko.computed(function () {
        return self.values.reduce(function(acc, val){
            return acc + (val.hasIteration() ? 2 : 1);
        },0);
    });
    self.getFields = function() {
        return [self];
    }

    self.serialized = ko.computed(function() {
        return [map(self.values, function (value) { return value.serialized(); })];
    });

    self.asHistory = ko.computed(function() {
        return [[self.label, "["+map(self.values, function (value) { return value.asHistory(); })+"]"]];
    });
}

function FilenameField(label, value) {
    var self = this;
    self.type = 3;
    self.label = label;
    self.value = value;
    self.url = value.replace(data_dir, "/dat");
    self.torender = value.replace(data_dir + "/", "");
    self.renderer = 'filename-row-template';
    self.getFields = function () { return [self]; }
    self.serialized = function () { return [self.value]; }
    self.asHistory = function() { return [[self.label, self.torender]]; }
}

function allequal(arr) {
    var a = arr[0];
    for (i = 1; i < arr.length; i++)
        if (a != arr[i])
            return false;
    return true;
}

function VectorCompactField(label, values) {
    var self = this;
    var isvector = allequal(values) == false;
    self.type = 2;
    self.label = label;
    self.renderer = 'vectorcompact-row-template';
    self.options = ["Constant", "Array"];
    self.isvector = ko.observable(isvector ? "Array" : "Constant");
    self.vector = $.map(values, function (value, i) {
        return new ScalarValue(value, label+'['+i+']');
    });
    self.scalar = self.vector[0];

    self.elements = ko.computed(function() {
        return self.isvector() == "Array" ? self.vector : [self.scalar];
    });
    self.spansize = ko.computed(function () {
        return 1 + self.elements().reduce(function(acc, val){
            return acc + (val.hasIteration() ? 2 : 1);
        },0);
    });

    self.at = function(i) {
        return self.elements()[i];
    }

    self.getFields = function() {
        return [self];
    }

    self.serialized = ko.computed(function() {
        if (self.isvector() == "Array")
            return [map(self.vector, function (value) { return value.serialized(); })];
        else {
            return [[self.scalar.serialized()]];
        }
    });
    self.asHistory = ko.computed(function() {
        var prefix = self.isvector() == "Array" ? "" : self.vector.length+"*";
        return [[self.label, prefix+"["+map(self.elements(), function (value) { return value.asHistory(); }) + "]"]];
    });
}


function iterkeys(arr) {
    acc = [];    
    for (i = 0; i < arr.length; i++)
        acc.push(arr[i][0]);
    return acc;
}

function lookup(arr, what) {
    for (i = 0; i < arr.length; i++)
        if (arr[i][0] == what)
            return arr[i][1];
    console.log('failed to find ' + what + ' in ' + arr);
    return undefined;
}

function flatten(arr){
    var res = arr.reduce(function(acc, val){
        return acc.concat(val.getFields());
    },[]);
    for (var i=0; i<res.length; i++)
        res[i].index = i;
    return res;
}

function serialized(arr){
    return arr.reduce(function(acc, val){
        return acc.concat(val.serialized());
    },[]);
}
function asHistory(arr){
    return arr.reduce(function(acc, val){
        return acc.concat(val.asHistory());
    },[]);
}

function EnumField(label, value, options_loaded, params, parent, entity) {
    // console.log(options)
    var self = this;

    var options = [];
    for (i = 0; i < options_loaded.length; i++) {
        var choice_label = options_loaded[i][0];
        var choice_params = choice_label==value ? params : options_loaded[i][1];
        options.push([choice_label, loadParams(choice_params, parent, entity.concat(label))]); // todo: pass true field name
    }

    self.type = options_loaded;
    self.label = label;
    self.value = ko.observable(value);
    self.options = iterkeys(options);
    self.renderer = 'enum-row-template';

    self.myParams = ko.computed(function() {
        return lookup(options, self.value());
    })

    self.getFields = ko.computed(function() {
        return [self].concat(flatten(self.myParams()));
    });

    self.serialized = ko.computed(function() {
        return [[self.value(), serialized(self.myParams())]];
    });

    self.asHistory = ko.computed(function() {
        var res = [[self.label, self.value()]];
        for (var i=0; i<self.myParams().length; i++) {
            var e = self.myParams()[i];
            res.concat(e.asHistory());
        }
        return res;
    });
}

function loadParams(raw, parent, entity) {
    return $.map(raw, function (e) {
        return (
            (e[1] == 0) ? new ScalarField(e[0], e[2], e[3], e[4], parent, entity, e[5]) :
            (e[1] == 1) ? new VectorField(e[0], e[2]) :
            (e[1] == 2) ? new VectorCompactField(e[0], e[2]) :
            (e[1] == 3) ? new FilenameField(e[0], e[2]) :
            new EnumField(e[0], e[2], e[1], e[3], parent, entity));
    });
}

function array_to_2d(src, len_1, len_2) {
    var res = []; 
    if (src.length != len_1*len_2)
        console.log("len_1="+len_1+"; len_2="+len_2+"; len(src)="+src.length);
    for (var i=0; i<len_1; i++) {
        res.push(src.splice(0, len_2));
    }
    return res;
}

function iterationRangOf(q) {
    var s = '"iterate"';
    var pos = q.indexOf(s, 0);
    if (pos == -1) return 0;
    pos = q.indexOf(s, pos+1);
    if (pos == -1) return 1;
    return 2;
}

function firstChild(e) {
    for (var j=0; j<e.childNodes.length; j++) {
        if (e.childNodes[j].nodeType == 1) {
            return e.childNodes[j];
        }
    }    
    return undefined;
}

function HistoryElement(root) {
    var self = this;
    self.myAsset = root.myAsset();
    self.myModel = root.myModel();
    self.myFamily = root.myFamilyId()[0];
    self.myOption = root.myOption()[0];
    self.myMethod = root.myMethod();
    self.myModelParams = asHistory(root.myModelParams());
    self.myOptionParams = asHistory(root.myOptionParams());
    self.myMethodParams = asHistory(root.myMethodParams());

    self.query = root.resultQuery();

    self.iterationRang = iterationRangOf(self.query);

    if (self.iterationRang == 0) {
        self.scalarResult = map(root.scalarResult(), function (e) { return [e[0](), e[1]]; });
    } else if (self.iterationRang == 1) {
        self.iteration1dLabel = root.iteration1dLabel; 
        self.iteration1dKinds = root.iteration1dKinds();

        self.iteration1dGraphsData = [];
        var kinds = self.iteration1dKinds;
        for (var k in kinds) {
            var kind = kinds[k];
            self.iteration1dGraphsData.push(kind);
        }
        self.renderGraph1d = root.renderGraph1d;
    } else if (self.iterationRang == 2) {
        self.iteration2dGraphsData = root.iteration2dGraphsData();
        self.renderGraph2d = root.renderGraph2d;
    }
}

function ModelView() {
    var self = this; 
    
    self.assets = get('assets');

    self.models = KsCachedMap(function(asset) {return 'models?a='+asset; });;
    self.families = KsCachedMap(function(model) {return 'families?m='+model; });;
    self.options = KsCachedMap(function(args){ return 'options?m='+args[0]+"&f="+args[1];});;
    self.methods = KsCachedMap(function(args){ return 'methods?m='+args[0]+"&f="+args[1]+"&o="+args[2];});

    self.model_html = KsCachedMap(function(m) {return 'model_html?m='+m; });;
    self.family_html = KsCachedMap(function(args) {return 'family_html?f='+args[0]+"&o="+args[1]; });;
    self.option_html = KsCachedMap(function(args){ return 'option_html?f='+args[0]+"&o="+args[1];});;
    self.method_html = KsCachedMap(function(args){ return 'method_html?m='+args[0]+"&f="+args[1]+"&meth="+args[2];});

    self.model_params = ParamCachedMap(self, 'model', function (model) { return 'model_params?m='+model; });
    self.option_params = ParamCachedMap(self, 'option', function (args) { return 'option_params?f='+args[0]+"&o="+args[1]; });
    self.method_params = ParamCachedMap(self, 'method', function (args) { return 'method_params?m='+args[0]+"&f="+args[1]+"&meth="+args[2]; });
    self.method_results = ResultCachedMap(function (args) { return 'method_results?m='+args[0]+"&f="+args[1]+"&meth="+args[2]; });

    //------------------------------------------------------- Asset

    self.myAsset = ko.observable(self.assets[0]);
    self.myAssetProxy = ko.computed({
        read: function () { return self.myAsset() },
        write: function (value) { 
            self.myAsset(value); 
        }
    });

    //-------------------------------------------------------- Model

    self.myModels = ko.computed(function(){
        return self.models.at(self.myAsset())();
    });
    
    self.myModel = ko.observable(self.myModels()[0]);

    self.myModelProxy = ko.computed({
        read: function () { return self.myModel() },
        write: function (value) { 
            self.myModel(value); 
        }
    });

    self.myModelHtml = ko.computed(function() {
        return self.model_html.at(self.myModel());
    });

    self.myModelParamsNF = ko.computed({
        read: function () { return self.model_params.at(self.myModel())(); },
        write: function (x) { self.model_params.set(self.myModel(), x); }
    });

    self.myModelParams = ko.computed(function () {
        return flatten(self.myModelParamsNF());
    });

    //---------------------------------------------------------- Family

    self.myFamilies = ko.computed(function(){
        var res = [self.families.at(self.myModel())(), self.myModel()];
        //console.log('computing myFamilies: myModel=' + self.myModel() + '; families <- ' + $.toJSON(res));
        return res;
    });

    self.myFamilyId = ko.observable([self.myFamilies()[0][0], self.myFamilies()[1]]);
    //self.myFamily = ko.observable(self.myFamilies()[0]);

    self.myFamilies.subscribe(function (families) {
        self.myFamilyId([families[0][0], families[1]]);
    })

    self.myFamilyProxy = ko.computed({
        read: function () { return self.myFamilyId()[0]; },
        write: function (value) { 
            //console.log("family="+value);
            //console.log("model="+self.myModel.peek());
            self.myFamilyId([value, self.myModel.peek()]); 
        }
    });

    //------------------------------------------------------------ Option

    self.myOptions = ko.computed(function(){
        var res = [self.options.at([self.myFamilyId()[1], self.myFamilyId()[0]])(), self.myFamilyId()[1]];
        //console.log('computing myOptions: myModel=' + self.myFamilyId()[1] + '; myFamily = ' + self.myFamilyId()[0] + '; families <- ' + $.toJSON(res));
        return res;
    });

    self.myOption = ko.observable([self.myOptions()[0][0], self.myFamilyId.peek()[1]]);

    self.myOptions.subscribe(function (options) {
        self.myOption([options[0][0],options[1]]);
    })

    self.myOptionProxy = ko.computed({
        read: function () { return self.myOption()[0]; },
        write: function (x) { self.myOption([x, self.myFamilyId.peek()[1]]); }
    })

    self.myOptionHtml = ko.computed(function() {
        return self.option_html.at([self.myFamilyId.peek()[0], self.myOption()[0]]);
    });

    self.myFamilyHtml = ko.computed(function() {
        return self.family_html.at([self.myFamilyId.peek()[0], self.myOption()[0]]);
    });
    self.myOptionParamsNF = ko.computed({
        read: function () { return self.option_params.at([self.myFamilyId.peek()[0], self.myOption()[0]])(); },
        write: function (x) { self.option_params.set([self.myFamilyId.peek()[0], self.myOption()[0]], x); }
    });
    self.myOptionParams = ko.computed(function () {
        return flatten(self.myOptionParamsNF());
    });

    //--------------------------------------------------------------- Method

    self.myMethods = ko.computed(function(){
        var methods = self.methods.at([self.myModel.peek(), self.myFamilyId.peek()[0], self.myOption()[0]])();
        //console.log('computing myMethods: myModel=' + self.myOption.peek()[1] +'myOption=' + self.myOption.peek()[0] + '; myFamily = ' + self.myFamilyId.peek()[0] + '; methods <- ' + $.toJSON(methods));
        return methods;
    });

    self.myMethod = ko.observable(self.myMethods()[0]);

    self.myMethodHtml = ko.computed(function() {
        return self.method_html.at([self.myFamilyId.peek()[1], self.myFamilyId.peek()[0], self.myMethod()]);
    });

    self.myMethodResults = ko.computed(function () {
        return self.method_results.at([self.myFamilyId.peek()[1], self.myFamilyId.peek()[0], self.myMethod()])();
    });

    self.myMethodParamsNF = ko.computed(function () {
        //console.log('method results updated: ' + self.myMethod());
        return self.method_params.at([self.myFamilyId.peek()[1], self.myFamilyId.peek()[0], self.myMethod()])();
    });

    self.myMethodParams = ko.computed(function () {
        return flatten(self.myMethodParamsNF());
    });
    //------------------------------------------------------------------- Result
    self.currentState = ko.computed(function(){
        return [
            [self.myAsset(), self.myModel(), serialized(self.myModelParamsNF())],
            [self.myFamilyId()[0], self.myOption()[0], serialized(self.myOptionParamsNF())],
            [self.myMethod(), serialized(self.myMethodParamsNF())]
        ];        
    })

    self.query = ko.computed(function () {
        return $.toJSON(self.currentState());
    });

    self.reload = function(path, value) {
        var query = $.toJSON([self.currentState(), [path, value]]);
        var response = get("adjust?"+query);
        console.log(response);
        var model = loadParams(response[0], self, ['model']);
        var option = loadParams(response[1], self, ['option']);
        var method = loadParams(response[2], self, ['method']);
        self.myModelParamsNF(model);
        self.myOptionParamsNF(option);
        self.myMethodParamsNF(method);
    }

    self.iterationRang = ko.computed(function() { 
        return iterationRangOf(self.query());
    })

    self.resultQuery = ko.observable("");

    self.resultIterationRang = ko.computed(function(){
        return iterationRangOf(self.resultQuery());
    })


    self.paramsAreOk = ko.computed(function() {
        var q = self.query();
        var s = '"_error_"';
        return q.indexOf(s, 0) == -1;        
    })

    self.history = ko.observableArray();
    self.history_element_to_push = undefined;

    self.scalarResult = ko.observable(undefined);
    self.errorResult = ko.observable(undefined);

    self.debug = ko.observable(false);
    self.graphSizeX = ko.observable(640);
    self.graphSizeY = ko.observable(384);
    self.inheritByDefault = ko.observable(true);
    self.trackSuffix = ko.observable(true);
    self.userSuffix = "";
    self.suffix = ko.computed({
        read: function () {
                return (self.trackSuffix() 
                            ? (lastLabel() && lastValue() 
                                    ? '(' + lastLabel() + '=' + lastValue() + ')' : '') 
                            : self.userSuffix);
        }, 
        write: function (x) {
            self.userSuffix = x;
            self.trackSuffix(false);
        }
    });

    //-------------------------------------------------------  Iteration 1D

    self.showGraph1d = ko.computed(function() {
        return self.resultIterationRang() == 1 && self.iteration1dKinds() != undefined;
    })

    self.iteration1dLabel = undefined;

    self.iteration1dKinds = ko.observable(undefined);

    self.iteration1dGraphsData = ko.computed(function() {
        var kinds = self.iteration1dKinds();        
        var graphs = [];
        if (kinds != undefined) {
            for (var k in kinds) {
                graphs.push(kinds[k]);
            }
        }

        return graphs;
    })

    self.renderGraph1d = function (elem, data) {
        //console.log("elem="+elem+",data="+data.label+", type="+map(elem, function (e) {return e.nodeType; }));
        for (var i=0; i<elem.length; i++) {
            var e = elem[i];
            if (e.nodeType==1) {
                var ee = firstChild(firstChild(firstChild(firstChild(e))));
                ee.style.width = self.graphSizeX()+'px';
                ee.style.height = self.graphSizeY()+'px';
                Flotr.draw(ee, data, {
                    legend : {
                        position : 'se',            // Position the legend 'south-east'.
                        backgroundColor : '#D2E8FF' // A light blue background color.
                    },
                    HtmlText : false,
                    spreadsheet: {
                      show: true,
                      tickFormatter: function(e) {
                          return e + "";
                      }
                    }

                });
            }
        }
    } 

    //-------------------------------------------------------- Iteration 2D

    self.showGraph2d = ko.computed(function() {
        return self.resultIterationRang() == 2;
    })

    self.iteration2dGraphsData = ko.observable(undefined);
    self.iteration2dLabels = undefined;

    self.renderGraph2d = function (elem, source) {
        //console.log("elem="+elem+",data="+data.label+", type="+map(elem, function (e) {return e.nodeType; }));
        for (var i=0; i<elem.length; i++) {
            var e = elem[i];
            if (e.nodeType==1) {
                var sizeX = parseInt(self.graphSizeX()) + 200;
                var sizeY = parseInt(self.graphSizeY()) + 200;

                e.style.width = sizeX+'px';
                e.style.height = sizeY+'px';
                
                var surfacePlot;
                var keys_1 = self.iteration2dLabels[0][1];
                var keys_2 = self.iteration2dLabels[1][1];
                var numRows = keys_1.length;
                var numCols = keys_2.length;

                // console.log("keys_1="+keys_1);
                // console.log("keys_2="+keys_2);

                var values = source.data;
                // for (var i=0; i<values.length; i++)
                //     console.log("values["+i+"]="+values[i]);

                var data = {nRows: numRows, nCols: numCols, formattedValues: values};

                surfacePlot = new SurfacePlot(e);

                // Don't fill polygons in IE. It's too slow.
                var fillPly = true;

                // Define a colour gradient.
                var colour1 = {red:0, green:0, blue:255};
                var colour2 = {red:0, green:255, blue:255};
                var colour3 = {red:0, green:255, blue:0};
                var colour4 = {red:255, green:255, blue:0};
                var colour5 = {red:255, green:0, blue:0};
                var colours = [colour1, colour2, colour3, colour4, colour5];

                // Axis labels.
                var xAxisHeader   = self.iteration2dLabels[0][0];
                var yAxisHeader   = self.iteration2dLabels[1][0];
                var zAxisHeader   = source.label;
                // console.log("xAxisHeader="+xAxisHeader);
                // console.log("yAxisHeader="+yAxisHeader);
                // console.log("zAxisHeader="+zAxisHeader);

                var renderDataPoints = false;
                var background = '#ffffff';
                var axisForeColour = '#000000';
                var hideFloorPolygons = true;
                var chartOrigin = {x: sizeX/2, y:sizeY/2};

                // Options for the basic canvas pliot.
                var basicPlotOptions = {fillPolygons: fillPly, tooltips: source.tooltips, renderPoints: renderDataPoints }

                // Options for the webGL plot.
                var xLabels = keys_1;
                var yLabels = keys_2;
                var glOptions = {xLabels: xLabels, yLabels: yLabels ,autoCalcZScale: true};

                // Options common to both types of plot.
                var options = {xPos: 0, yPos: 0, width: sizeX, height: sizeY, colourGradient: colours,
                 xTitle: xAxisHeader, yTitle: yAxisHeader, zTitle: zAxisHeader,
                 backColour: background, axisTextColour: axisForeColour, hideFlatMinPolygons: hideFloorPolygons, origin: chartOrigin};

                surfacePlot.draw(data, options, basicPlotOptions, glOptions);
            }
        }
    } 

    self.hasResult = ko.computed(function() {
        return self.errorResult() || self.scalarResult() || self.iteration1dKinds() || self.iteration2dGraphsData();
    });

    self.computing = ko.computed(function () { 
        return self.query() == self.resultQuery() && !self.hasResult();
    });

    var historyPusher = ko.computed(function () {
        if (self.query() != self.resultQuery()) {
            if (self.hasResult()) {
                if (self.history_element_to_push) {
                    self.history.unshift(self.history_element_to_push);
                    self.history_element_to_push = undefined;
                }
                self.clearResult();
                self.resultQuery("");
            }
        }
        return 0;
    })

    self.clearResult = function() {
        self.scalarResult(undefined);
        self.errorResult(undefined);
        self.iteration1dLabel = undefined;
        self.iteration1dKinds(undefined);
        self.iteration2dGraphsData(undefined);
        self.iteration2dLabels = undefined;
    }

    self.updateResult = function(result) {
        if (result[0][0] == "Exception") {
            self.errorResult(result[0][1][0]);
        } else {

            if (self.resultIterationRang() == 0) {
                self.scalarResult(map(result, function (val, i) { return [self.myMethodResults()[i].label, val[1]]; })); 
            } else if (self.resultIterationRang() == 1) {
                var kinds = {};
                var src = result;
                self.iteration1dLabel = result[0][0];
                for (var i = 1; i<src.length; i++) {
                    var meta = self.myMethodResults()[i-1];
                    if (meta.visible()) {
                        var s = src[i];            
                        if (typeof(s[1][0]) == "number") {
                            var d = {
                                data: map(s[1], function(x,j) { return [src[0][1][j], x]; }),
                                label: meta.label() + self.suffix(),
                                inherit: ko.observable(self.inheritByDefault())
                            }
                            if (kinds[meta.kind] == undefined)
                                kinds[meta.kind] = [];
                            kinds[meta.kind].push(d);
                        } else {
                            for (var ii=0; ii < s[1][0].length; ii++) {
                                var d = {
                                    data: map(s[1], function(x,j) { return [src[0][1][j], x[ii]]; }),
                                    label: meta.label()+'['+ii+']'+ self.suffix(),
                                    inherit: ko.observable(self.inheritByDefault())
                                }
                                if (kinds[meta.kind] == undefined)
                                    kinds[meta.kind] = [];
                                kinds[meta.kind].push(d);
                            }
                        }
                    }
                }
                var history = self.history();            
                for (var h_idx=0; h_idx<history.length; h_idx++ ) {
                    var h = history[h_idx];
                    if (h.iteration1dLabel == self.iteration1dLabel) {
                        var hKinds = h.iteration1dKinds; 
                        for (var k_idx in hKinds) {
                            if (kinds[k_idx]==undefined) 
                                kinds[k_idx] = [];
                            var hKind = hKinds[k_idx];
                            for (var i=0; i<hKind.length; i++) {
                                if (hKind[i].inherit()) {
                                    var duplicata = false;
                                    for (var z=0; z<kinds[k_idx].length; z++) {
                                        if (kinds[k_idx][z] == hKind[i]) {
                                            duplicata = true;
                                            break;
                                        }
                                    }
                                    if (!duplicata)
                                        kinds[k_idx].push(hKind[i]);        
                                }
                            }
                        }           
                    } 
                }
                self.iteration1dKinds(kinds);
            } else if (self.resultIterationRang() == 2) {
                var graphs = [];
                var src = result;
                self.iteration2dLabels = [result[0], result[1]];
                for (var i = 2; i<src.length; i++) {
                    console.log("src="+src);         
                    var len_1 = src[0][1].length;
                    var len_2 = src[1][1].length;
                    var meta = self.myMethodResults()[i-2];
                    if (meta.visible()) {
                        var s = src[i];   
                        console.log("s="+s);         
                        if (typeof(s[1][0]) == "number") {
                            var d = {
                                tooltips: map(s[1], function(x) { return ""+x; }),
                                data: array_to_2d(s[1], len_1, len_2),
                                label: meta.label()
                            }
                            graphs.push(d);
                        } else {
                            for (var ii=0; ii < s[1][0].length; ii++) {
                                var d = {
                                    tooltips: map(s[1], function(x,j) { return ""+x[ii]; }),
                                    data: array_to_2d(map(s[1], function(x,j) { return x[ii]; }), len_1, len_2),
                                    label: meta.label()+'['+ii+']'
                                }
                                graphs.push(d);
                            }
                        }
                    }
                }
                self.iteration2dGraphsData(graphs);
            }
            self.history_element_to_push = new HistoryElement(self);            
        }
    }    
}

mv = new ModelView();

ko.applyBindings(mv);

$('#Compute').click(function() { 
    var my_query = mv.query();
    mv.clearResult();
    mv.resultQuery(my_query);
    $.getJSON('api.ks/compute?'+my_query, function (data) {
        console.log(data);
        if (my_query == mv.query()) {
            mv.updateResult(data);
        }
    }); 
});

$('#ClearHistory').click(function() {
    mv.history.removeAll();
});