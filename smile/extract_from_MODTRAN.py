# Author: Unknown
# This file extracts an array of data from MODTRAN corresponding to a specific feature

def extract_from_MODTRAN(json: dict, feature_wavelength: int, deeper_by = 0) -> np.array:
  """This function returns an array of data points corresposnding to the specific feature.
  
  Args:
      json = the docs_json variable copied from the website source
      feature_wavelength = the feature wavelength around which data is to be extracted
      deeper_by = the amount by which the y-axis is to be lowered in order to fit in with the satellite data
      
  Outputs:
      np.array consisting of x and y axis data points, wavelength in nm on x-axis and transmittance on the y-axis"""
  key = list(json.keys())[0]
  ind = 0
  for i in range(len(json[key]['roots']['references'])):
      if 'data' in json[key]['roots']['references'][i]['attributes']:
          ind = i

  xpoint = np.array(
      json[key]['roots']['references'][ind]['attributes']['data']['x'])
  ypoints = np.array(
      json[key]['roots']['references'][ind]['attributes']['data']['y'])
  xpoints = np.array([d * 1000 for d in xpoint])

  n1, n2 = feature_wavelength - 50, feature_wavelength + 50
  xp = np.array([d for d in xpoints if n1 < d < n2])
  pos = [list(xpoints).index(d) for d in list(xpoints) if n1 < d < n2]
  optimised = []
  for i in pos:
      if ypoints[i] < 0.87:
          optimised.append(ypoints[i] - deeper_by)
      else:
          optimised.append(ypoints[i])
  yp = np.array(optimised)
  plt.plot(xp, yp)
  plt.show()
  
  return np.array([xp, yp])


doc = {"56b32a2b-fe9c-4115-8c3f-c692609fd3e1":{"roots":{"references":[{"attributes":{"text":""},"id":"51570","type":"Title"},
    {"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":
        {"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"render_mode":"css","right_units":"screen","top_units":"screen"}
        ,"id":"51577","type":"BoxAnnotation"},{"attributes":{},"id":"51550","type":"WheelZoomTool"},{"attributes":{},"id":"51575","type":"Selection"},{"attributes":
        {"axis_label":"WAVELENGTH (microns)","formatter":{"id":"51572","type":"BasicTickFormatter"},"ticker":{"id":"51540","type":"BasicTicker"}},"id":"51539","type":"LinearAxis"},
        {"attributes":{},"id":"51537","type":"LinearScale"},{"attributes":{},"id":"51535","type":"LinearScale"},{"attributes":{},"id":"51552","type":"ResetTool"},{"attributes":
        {"callback":'null'},"id":"51533","type":"DataRange1d"},{"attributes":{"callback":'null'},"id":"51531","type":"DataRange1d"},{"attributes":{"below":[{"id":"51539","type":"LinearAxis"}]
        ,"center":[{"id":"51543","type":"Grid"},{"id":"51548","type":"Grid"}],"left":[{"id":"51544","type":"LinearAxis"}],"plot_height":400,"plot_width":400,"renderers":[{"id":"51567","type":"GlyphRenderer"}]
        ,"title":{"id":"51570","type":"Title"},"toolbar":{"id":"51556","type":"Toolbar"},"x_range":{"id":"51531","type":"DataRange1d"},"x_scale":{"id":"51535","type":"LinearScale"},"y_range":{"id":"51533","type":"DataRange1d"}
        ,"y_scale":{"id":"51537","type":"LinearScale"}},"id":"51530","subtype":"Figure","type":"Plot"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"}
        ,"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"render_mode":"css","right_units":"screen","top_units":"screen"}
        ,"id":"51578","type":"BoxAnnotation"},{"attributes":{"callback":'null'},"id":"51555","type":"HoverTool"},{"attributes":{"callback":'null',"overlay":{"id":"51578","type":"BoxAnnotation"}}
        ,"id":"51554","type":"BoxSelectTool"},{"attributes":{"active_drag":"auto","active_inspect":"auto","active_multi":'null',"active_scroll":"auto","active_tap":"auto","tools":[{"id":"51549","type":"PanTool"}
        ,{"id":"51550","type":"WheelZoomTool"},{"id":"51551","type":"BoxZoomTool"},{"id":"51552","type":"ResetTool"},{"id":"51553","type":"SaveTool"},{"id":"51554","type":"BoxSelectTool"}
        ,{"id":"51555","type":"HoverTool"}]},"id":"51556","type":"Toolbar"},{"attributes":{"overlay":{"id":"51577","type":"BoxAnnotation"}},"id":"51551","type":"BoxZoomTool"},{"attributes":{}
        ,"id":"51576","type":"UnionRenderers"},{"attributes":{},"id":"51553","type":"SaveTool"},{"attributes":{},"id":"51574","type":"BasicTickFormatter"},{"attributes":{"callback":'null'
        ,"data":{"x":[1.1,1.101555,1.10311,1.104665,1.10622,1.107775,1.10933,1.110885,1.11244,1.113995,1.11555,1.117105,1.11866,1.120215,1.12177,1.123325,1.12488,1.126435,1.12799,1.129545,1.1311,1.132655,1.13421,1.135765,1.13732,1.138875,1.14043,1.141985,1.14354,1.145095,1.14665,1.148205,1.14976,1.151315,1.15287,1.154425,1.15598,1.157535,1.15909,1.160645,1.1622,1.163755,1.16531,1.166865,1.16842,1.169975,1.17153,1.173085,1.17464,1.176195,1.17775,1.179305,1.18086,1.182415,1.18397,1.185525,1.18708,1.188635,1.19019,1.191745,1.1933,1.194855,1.19641,1.197965,1.19952,1.201075,1.20263,1.204185,1.20574,1.207295,1.20885,1.210405,1.21196,1.213515,1.21507,1.216625,1.21818,1.219735,1.22129,1.222845,1.2244,1.225955,1.22751,1.229065,1.23062,1.232175,1.23373,1.235285,1.23684,1.238395,1.23995,1.241505,1.24306,1.244615,1.24617,1.247725,1.24928,1.250835,1.25239,1.253945,1.2555,1.257055,1.25861,1.260165,1.26172,1.263275,1.26483,1.266385,1.26794,1.269495,1.27105,1.272605,1.27416,1.275715,1.27727,1.278825,1.28038,1.281935,1.28349,1.285045,1.2866,1.288155,1.28971,1.291265,1.29282,1.294375,1.29593,1.297485,1.29904,1.300595,1.30215,1.303705,1.30526,1.306815,1.30837,1.309925,1.31148,1.313035,1.31459,1.316145,1.3177,1.319255,1.32081,1.322365,1.32392,1.325475,1.32703,1.328585,1.33014,1.331695,1.33325,1.334805,1.33636,1.337915,1.33947,1.341025,1.34258,1.344135,1.34569,1.347245,1.3488,1.350355,1.35191,1.353465,1.35502,1.356575,1.35813,1.359685,1.36124,1.362795,1.36435,1.365905,1.36746,1.369015,1.37057,1.372125,1.37368,1.375235,1.37679,1.378345,1.3799,1.381455,1.38301,1.384565,1.38612,1.387675,1.38923,1.390785,1.39234,1.393895,1.39545,1.397005,1.39856,1.400115,1.40167,1.403225,1.40478,1.406335,1.40789,1.409445,1.411,1.412555,1.41411,1.415665,1.41722,1.418775,1.42033,1.421885,1.42344,1.424995,1.42655,1.428105,1.42966,1.431215,1.43277,1.434325,1.43588,1.437435,1.43899,1.440545,1.4421,1.443655,1.44521,1.446765,1.44832,1.449875,1.45143,1.452985,1.45454,1.456095,1.45765,1.459205,1.46076,1.462315,1.46387,1.465425,1.46698,1.468535,1.47009,1.471645,1.4732,1.474755,1.47631,1.477865,1.47942,1.480975,1.48253,1.484085,1.48564,1.487195,1.48875,1.490305,1.49186,1.493415,1.49497,1.496525,1.49808,1.499635,1.50119,1.502745,1.5043,1.505855,1.50741,1.508965,1.51052,1.512075,1.51363,1.515185,1.51674,1.518295,1.51985,1.521405,1.52296,1.524515,1.52607,1.527625,1.52918,1.530735,1.53229,1.533845,1.5354,1.536955,1.53851,1.540065,1.54162,1.543175,1.54473,1.546285,1.54784,1.549395,1.55095,1.552505,1.55406,1.555615,1.55717,1.558725,1.56028,1.561835,1.56339,1.564945,1.5665,1.568055,1.56961,1.571165,1.57272,1.574275,1.57583,1.577385,1.57894,1.580495,1.58205,1.583605,1.58516,1.586715,1.58827,1.589825,1.59138,1.592935,1.59449,1.596045,1.5976,1.599155,1.60071,1.602265,1.60382,1.605375,1.60693,1.608485,1.61004,1.611595,1.61315,1.614705,1.61626,1.617815,1.61937,1.620925,1.62248,1.624035,1.62559,1.627145,1.6287,1.630255,1.63181,1.633365,1.63492,1.636475,1.63803,1.639585,1.64114,1.642695,1.64425,1.645805,1.64736,1.648915,1.65047,1.652025,1.65358,1.655135,1.65669,1.658245,1.6598,1.661355,1.66291,1.664465,1.66602,1.667575,1.66913,1.670685,1.67224,1.673795,1.67535,1.676905,1.67846,1.680015,1.68157,1.683125,1.68468,1.686235,1.68779,1.689345,1.6909,1.692455,1.69401,1.695565,1.69712,1.698675,1.70023,1.701785,1.70334,1.704895,1.70645,1.708005,1.70956,1.711115,1.71267,1.714225,1.71578,1.717335,1.71889,1.720445,1.722,1.723555,1.72511,1.726665,1.72822,1.729775,1.73133,1.732885,1.73444,1.735995,1.73755,1.739105,1.74066,1.742215,1.74377,1.745325,1.74688,1.748435,1.74999,1.751545,1.7531,1.754655,1.75621,1.757765,1.75932,1.760875,1.76243,1.763985,1.76554,1.767095,1.76865,1.770205,1.77176,1.773315,1.77487,1.776425,1.77798,1.779535,1.78109,1.782645,1.7842,1.785755,1.78731,1.788865,1.79042,1.791975,1.79353,1.795085,1.79664,1.798195,1.79975,1.801305,1.80286,1.804415,1.80597,1.807525,1.80908,1.810635,1.81219,1.813745,1.8153,1.816855,1.81841,1.819965,1.82152,1.823075,1.82463,1.826185,1.82774,1.829295,1.83085,1.832405,1.83396,1.835515,1.83707,1.838625,1.84018,1.841735,1.84329,1.844845,1.8464,1.847955,1.84951,1.851065,1.85262,1.854175,1.85573,1.857285,1.85884,1.860395,1.86195,1.863505,1.86506,1.866615,1.86817,1.869725,1.87128,1.872835,1.87439,1.875945,1.8775,1.879055,1.88061,1.882165,1.88372,1.885275,1.88683,1.888385,1.88994,1.891495,1.89305,1.894605,1.89616,1.897715,1.89927]
        ,"y":[0.8792122,0.87951809,0.87979984,0.8800807,0.88036197,0.88064206,0.88092035,0.88119292,0.88145751,0.88172334,0.88199842,0.88227504,0.8825528,0.88282955,0.88310045,0.88336861,0.883636,0.88390303,0.88416421,0.88441479,0.88466829,0.88493192,0.88519603,0.88546419,0.88573396,0.88599807,0.88625753,0.88651377,0.88676804,0.88702208,0.88727576,0.88752824,0.88777858,0.88802856,0.8882792,0.88852775,0.8887772,0.88902688,0.88927418,0.88951904,0.88976216,0.89000553,0.89024889,0.8904919,0.89073354,0.89097363,0.89121205,0.89144933,0.89168578,0.89189613,0.89206636,0.8922416,0.89244264,0.89266735,0.89290452,0.89311963,0.89331692,0.89351183,0.89371282,0.8939448,0.89419818,0.89443237,0.89464426,0.89484656,0.89503217,0.89437789,0.88965058,0.88491565,0.8861149,0.88809866,0.88949418,0.89220738,0.89431441,0.89536625,0.89576858,0.89397675,0.89009511,0.88903511,0.89054483,0.89117754,0.89208919,0.89341229,0.89405733,0.89430881,0.89426869,0.89401424,0.89368516,0.89334339,0.89298552,0.89254636,0.89203763,0.89143455,0.89069831,0.88978678,0.88858646,0.88696843,0.88481104,0.88200504,0.87846845,0.87418193,0.8692736,0.86300159,0.85310549,0.83714062,0.81554389,0.79400492,0.78410172,0.77225709,0.73569942,0.7333672,0.79001898,0.82858139,0.83859557,0.84971184,0.86103398,0.8717342,0.88055569,0.88724375,0.89181107,0.89479512,0.89713383,0.89907897,0.90041876,0.90152031,0.90269238,0.90360624,0.9043951,0.90522355,0.90591681,0.90645963,0.90692753,0.90733469,0.90770501,0.9080686,0.90841025,0.90868104,0.90866834,0.90804881,0.90755081,0.9080016,0.90851128,0.90868986,0.90921491,0.90997565,0.91059303,0.91103768,0.91133875,0.91154522,0.91172069,0.91188598,0.91202396,0.91203219,0.91184235,0.91176265,0.91204304,0.91231728,0.91242844,0.91267741,0.91306019,0.91340011,0.91365093,0.9138329,0.91400051,0.91417587,0.91434926,0.91451627,0.91467738,0.91483587,0.91499519,0.91515541,0.91531217,0.91546386,0.91560763,0.91574389,0.91589284,0.91605926,0.91621655,0.91636515,0.91652822,0.9167006,0.91686481,0.91701841,0.91716814,0.91732031,0.91747302,0.91762424,0.91777581,0.91792732,0.91807741,0.91822755,0.91837674,0.91852266,0.91866875,0.91881794,0.91896671,0.91911268,0.91925895,0.919406,0.91955298,0.91969848,0.9198423,0.91998559,0.92012858,0.92027187,0.9204151,0.92055732,0.92069852,0.92083871,0.92097801,0.92111552,0.9212451,0.91820306,0.85574394,0.74894643,0.74247146,0.80480218,0.81727791,0.81055391,0.82130951,0.84655082,0.87455076,0.8914997,0.90153641,0.90921509,0.9148612,0.91867208,0.92084777,0.92205203,0.9228133,0.92334312,0.92370635,0.92396188,0.92415851,0.92432648,0.92447782,0.92461741,0.92475098,0.92481869,0.92378521,0.9218145,0.92172211,0.92317683,0.92357886,0.92339861,0.92373735,0.92448193,0.92520696,0.92566025,0.92594701,0.92618257,0.92638814,0.92651719,0.92648327,0.92637479,0.92633927,0.926287,0.92642921,0.92693663,0.92736846,0.9275772,0.9277364,0.92787528,0.92800266,0.92812353,0.92824221,0.92835951,0.9284569,0.92847592,0.92848867,0.92864043,0.92882562,0.92889982,0.92882615,0.92855382,0.92829108,0.92827237,0.92811114,0.9269321,0.92406058,0.92080468,0.92072189,0.92369473,0.92468601,0.92323232,0.92350107,0.92598635,0.92852718,0.93012494,0.930933,0.93133742,0.93155801,0.93172395,0.93188095,0.93202835,0.93215656,0.93227708,0.93241316,0.93257082,0.9326883,0.93242162,0.92985982,0.91816264,0.89087015,0.8621009,0.85911453,0.88017362,0.89099431,0.8812924,0.877738,0.88992149,0.90577054,0.91881865,0.9274832,0.93194783,0.93384212,0.9346205,0.93497062,0.9351247,0.93503469,0.9336381,0.92589563,0.90581155,0.87835842,0.8619777,0.87336707,0.89425963,0.89229488,0.87980533,0.88150573,0.89506561,0.91044974,0.92159277,0.92785299,0.93152881,0.93398345,0.93558627,0.93657327,0.93716955,0.93745601,0.93747365,0.93744278,0.93760872,0.93792379,0.93813008,0.93806058,0.93755722,0.93617886,0.93375492,0.93131071,0.93079758,0.93315268,0.93533236,0.93427753,0.93271053,0.93352121,0.93564773,0.93756115,0.93877745,0.93942058,0.93985665,0.94027954,0.94049317,0.94044155,0.94044274,0.9406445,0.94089633,0.94107777,0.94123888,0.94142985,0.94159418,0.94172019,0.94186413,0.94202119,0.94216037,0.94229311,0.94243294,0.94257903,0.94272149,0.94285482,0.94297808,0.94308543,0.94317538,0.94326442,0.94337445,0.9435041,0.94362307,0.94370532,0.94376153,0.94383699,0.94396591,0.94411832,0.94423956,0.94433552,0.94444239,0.9445591,0.94467515,0.94479585,0.94491804,0.94503373,0.94514358,0.9452495,0.94535387,0.9454577,0.94556153,0.94566417,0.94576567,0.94586718,0.94596928,0.9460721,0.94617379,0.94627392,0.94637424,0.94647527,0.94657648,0.94667715,0.94677752,0.946877,0.94697648,0.94707507,0.94715971,0.94721019,0.94726688,0.94738328,0.94751674,0.94761503,0.94770002,0.94779348,0.94785428,0.94786823,0.9479183,0.94805902,0.948228,0.94833881,0.94840646,0.94849092,0.94859743,0.94871306,0.94882792,0.94893652,0.94904238,0.94912392,0.94917601,0.94925445,0.94938046,0.94950539,0.94959968,0.94968015,0.94976813,0.94986385,0.94996083,0.95005751,0.95015103,0.9502393,0.95032519,0.95041591,0.95051146,0.95060492,0.95069331,0.95077986,0.95086741,0.95095575,0.9510448,0.95113409,0.95122278,0.9513101,0.95139706,0.95148343,0.9515698,0.95165575,0.95174056,0.95182502,0.95190954,0.95199376,0.95207804,0.95216262,0.9522469,0.95233077,0.95241398,0.95249623,0.95257872,0.95266134,0.95274383,0.95282573,0.95290715,0.95298839,0.95306945,0.95315075,0.95323205,0.95331341,0.95339429,0.95347422,0.95355022,0.95362258,0.95369524,0.95372474,0.95370543,0.95345885,0.9522835,0.9500466,0.94784063,0.94735628,0.94652188,0.93950015,0.93491262,0.94197011,0.94910157,0.95087361,0.95159394,0.95201236,0.95276535,0.95369548,0.95426565,0.95457137,0.95476937,0.95489198]}
        ,"selected":{"id":"51575","type":"Selection"},"selection_policy":{"id":"51576","type":"UnionRenderers"}},"id":"51564","type":"ColumnDataSource"},{"attributes":{"line_color":"#1f77b4","x":{"field":"x"},"y":{"field":"y"}},"id":"51565","type":"Line"}
        ,{"attributes":{"line_alpha":0.1,"line_color":"#1f77b4","x":{"field":"x"},"y":{"field":"y"}},"id":"51566","type":"Line"},{"attributes":{"data_source":{"id":"51564","type":"ColumnDataSource"}
        ,"glyph":{"id":"51565","type":"Line"},"hover_glyph":'null',"muted_glyph":'null',"nonselection_glyph":{"id":"51566","type":"Line"},"selection_glyph":'null',"view":{"id":"51568","type":"CDSView"}}
        ,"id":"51567","type":"GlyphRenderer"},{"attributes":{"source":{"id":"51564","type":"ColumnDataSource"}},"id":"51568","type":"CDSView"},{"attributes":{"dimension":1,"ticker":{"id":"51545","type":"BasicTicker"}}
        ,"id":"51548","type":"Grid"},{"attributes":{},"id":"51549","type":"PanTool"},{"attributes":{"axis_label":"TRANSMITTANCE","formatter":{"id":"51574","type":"BasicTickFormatter"},"ticker":{"id":"51545","type":"BasicTicker"}}
        ,"id":"51544","type":"LinearAxis"},{"attributes":{},"id":"51545","type":"BasicTicker"},{"attributes":{"ticker":{"id":"51540","type":"BasicTicker"}},"id":"51543","type":"Grid"},{"attributes":{},"id":"51540","type":"BasicTicker"}
        ,{"attributes":{},"id":"51572","type":"BasicTickFormatter"}],"root_ids":["51530"]},"title":"Bokeh Application","version":"1.4.0"}}


extract_from_MODTRAN(doc, 1450)