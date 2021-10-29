Models = {
    # Single componet models (note that "step": a cutoff powerlaw and "grbm": the Band function)
    "Power law (PL)": [['powerlaw'], 
                       {1: '-1.5'}], 
    "Blackbody (BB)": [['bbody'], 
                       {1: '100, 0.001, 10, 10, 1000, 1000', 2: '1'}], 
    "Cutoff PL (CPL)": [['cutoffpl'], 
                        {1: '-1, 0.001, -5, -10, 5, 10', 2: '100,0.1,10,10,100000,100000'}],
    "Cutoff PL (2pts approx.)": [['step'], 
                                 {1: '1, 0.001, 0, 0, 2, 10', 2: '200,0.01,10,10,1000000,1000000', 3:'1, 0.01,0.01,0.001, 10, 10'}],
    "Band": [['grbm'], 
             {1: '-0.6, 0.001, -5, -10, 2, 5', 2: '-2.3', 3: '500,0.1,10,10,10000,10000'}], 
    "Band with cutoff (BandwC)": [['grbm*highecut'], 
                         {1: '-0.7', 2: '-2.3', 3: '300,0.1,10,10,10000,10000', 
                          5: '10000, 1, 1000, 1000, 1000000, 1000000', 6: '10000, 1, 10, 10, 1000000, 1000000'}], 
    "Broken PL (BknPL)": [['bknpower'],
              {1: '-1', 2:'100000, 0.001, 100, 100, 10000000, 1000000', 3: '-2.2'}],
    "BknPL with cutoff (BknPLwC)": [['bknpower*highecut'], 
                          {1: '-1.0', 2: '300,0.1,10,10,10000,10000', 3: '-2.3', 
                           5: '10000, 1, 1000, 1000, 1000000, 1000000', 6: '10000, 1, 10, 10, 1000000, 1000000'}],
    "Two BknPL (Bkn2PL)": [['bkn2power'],
                  {1: '-0.7', 2: '-1.5', 3:'-2.5', 4: '300,0.1,10,10,10000,10000', 5:'3000, 1, 1000, 1000, 1000000, 1000000'}],

    # Two component models
    "PL + BB": [['powerlaw', 'bbody'], 
                {1: '-1.5', 
                 3: '40,0.1,10,10,10000,100000', 4: '1'}], 
    "CPL + PL": [['cutoffpl', 'powerlaw'],
                 {1: '-0.7, 0.001, -3, -5, 3, 5', 2: '1500,1,10,10,100000,100000', 3:'1, 0.001, 1e-06, 1e-08, 100, 100',
                  4: '-1.7', 5:'1, 0.001, 0, 0'}], 
    "CPL2 + PL": [['step', 'powerlaw'],
                  {1: '0.7, 0.001, 0, 0, 2, 5', 2: '300,0.01,10,10,1000000,1000000', 3:'1, 0.01,0.01,0.001, 10, 10',
                   4: '-1.8', 5:'1, 0.001, 0, 0'}], 
    "CPL + BB": [['cutoffpl', 'bbody'],
                 {1: '-1.5, 0.001, -5, -10, 5, 10', 2: '2000,0.01,10,10,1000000,1000000', 
                  4: '40,0.01,10,10,1000,1000', 5: '1'}], 
    "Two CPLs": [['step', 'cutoffpl'], 
                 {1: '1.5, 0.001, 0, 0, 5, 10', 2: '20000,0.01,10,10,1000000,1000000', 3:'1, 0.01,0.01,0.001, 10, 10', 
                  4: '-0.5, 0.001, -3, -10, 3, 10', 5: '200,0.01,10,10,1000000,1000000'}],
    "Band + PL": [['grbm', 'powerlaw'],
                  {1: '-0.7, 0.01', 2: '-2.3, 0.01,-10,-10,-1,-1', 3: '300,0.1,10,10,100000,100000', 
                   5: '-1.5, 0.01, -3, -5, 3, 5'}],
    "Band + BB": [['grbm', 'bbody'], 
                  {1: '-1', 2: '-2.2', 3: '1000,0.1,100,100,100000,100000', 
                   5: '20,0.001,20,20,1000,1000', 6: '1'}], 
    "Band + CPL": [['grbm', 'cutoffpl'],
                   {1: '-1, 0.01, -5, -10, 5, 10', 2: '-2.4, 0.01,-10,-10,-1,-1', 3: '500,0.1,10,10,1000000,10000000', 
                    5: '-0.8, 0.001, -5, -10, 5, 10', 6: '2000,0.1,10,10,2000000,2000000'}],
    "Band + CPL2": [['grbm', 'step'],
                    {1: '-1, 0.01, -5, -10, 5, 10', 2: '-2.4, 0.01,-10,-10,-1,-1', 3: '500,0.1,10,10,1000000,10000000',
                     5: '1.5, 0.001, 0, 0, 5, 10', 6: '200000,0.01,10,10,1000000,1000000', 7:'1, 0.001,0.001,0.001, 10, 10',}],
    "Band + CPL2": [['grbm', 'bknpower'],          
                    {1: '-1, 0.01, -5, -10, 5, 10', 2: '-2.4, 0.01,-10,-10,-1,-1', 3: '200,0.1,10,10,1000000,10000000', 
                     5: '1, 0.001, -5, -10, 15, 20', 6: '2000,0.1,100,100,5000000,5000000', 7:'-2.2'}],
    "BandwC + PL": [['grbm*highecut', 'powerlaw'], 
                    {1: '-1.0', 2: '-2.3', 3: '300,0.1,10,10,10000,10000', 
                     5:'10000, 1, 1000, 1000, 1000000, 1000000', 6:'100000, 1, 10, 10, 1000000, 1000000', 
                     7: '-2.0,0.01,-3,-5,3,5', 8: '1'}],
    "BandwC + BB": [['grbm*highecut', 'bbody'],    
                    {1: '-1.0', 2: '-2.3', 3: '300,0.1,10,10,10000,10000', 
                     5:'10000, 1, 1000, 1000, 1000000, 1000000', 6:'100000, 1, 10, 10, 1000000, 1000000', 
                     7: '300,0.001,20,20,1000,1000', 8: '1'}], 
    "BknPL + PL": [['bknpower', 'powerlaw'],      
                   {1: '-0.7, 0.01', 2: '-2.3, 0.01,-10,-10,-1,-1', 3: '300,0.1,10,10,100000,100000',
                    5: '-1.5, 0.01, -3, -5, 3, 5'}],
    "BknPL + BB": [['bknpower', 'bbody'],
                   {1: '-1', 2: '350,0.1,100,100,100000,100000', 3: '-2.2', 
                    5: '300,0.001,20,20,1000,1000', 6: '1'}], 
    "BknPL + CPL": [['bknpower', 'cutoffpl'],      
                    {1: '-0.7, 0.001, -10, -10, 20, 20', 2:'-2.2', 3:'30000, 0.1, 100, 100, 100000, 100000', 4: '-1, 0.001, -3, -5, 3, 5', 
                     5: '200,1,100,100,100000,100000', 6:'1, 0.001, 1e-06, 1e-08, 100, 100'}], 
    "BknPL + CPL2": [['bknpower', 'step'],
                     {1: '-1, 0.01, -5, -10, 5, 10', 2: '500,0.1,10,10,1000000,10000000', 3: '-2.4, 0.01,-10,-10,-1,-1', 4: '1, 0.001, -2, -10, 2, 10', 
                      5: '20000,0.01,10,10,1000000,1000000', 6:'1, 0.01,0.01,0.001, 10, 10'}],
    "Bkn2PL + PL": [['bkn2power', "powerlaw"],     
                    {1: '-0.7', 2: '-1.5', 3:'-2.5', 4: '300,0.1,10,10,10000,10000', 5:'3000, 1, 100, 100, 1000000, 1000000', 
                     7:"-1.8"}], 
    
    # Three component models
    "CPL + PL + BB": [['cutoffpl', 'powerlaw', 'bbody'],  
                      {1: '-0.7', 2: '2300,1,1,10, 100000,100000', 3:'0.179, 0.001, 0, 0, 100, 100', 
                       4: '-1.5', 5:'11.4, 0.001, 0, 0', 
                       6: '350,0.1, 1, 10, 10000,10000', 7:'17.2, 0.001, 0, 0'}], 
    "Band + PL + BB": [['grbm', 'powerlaw', 'bbody'],
                       {1: '-1', 2: '-2.2', 3: '300,1,1,10, 100000,100000', 4:'1, 0.001, 0, 0, 100, 100', 
                        5: '-1.5', 6:'1, 0.001, 0, 0', 
                        7: '100,0.1, 1, 10, 10000,10000', 8:'1, 0.001, 0, 0'}], 
    "CPL + CPL + PL": [['step', 'cutoffpl', 'powerlaw'],   
                       {1: '0.7', 2: '3000,0.01,10,10,1000000,1000000', 3:'1,0.01,0.01,0.001, 10, 10', 
                        4:'0.5', 5: '2000,1,1,10, 100000,100000', 
                        7: '-1.5'}],
    "CPL + CPL + BB": [['step', 'cutoffpl', 'bbody'], 
                       {1: '1.6, 0.001, 0, 0, 5, 10', 2: '2000,0.1,10,10,1000000,1000000', 3: '1, 0.01,0.01,0.001, 10, 10', 
                        4: '-0.5, 0.001, -3, -10, 3, 10', 5: '500,0.01,10,10,1000000,1000000', 
                        7: '100,0.1, 1, 10, 10000,10000', 8:'1, 0.001, 0, 0'}],
    "Band + CPL + BB": [['grbm', 'cutoffpl', 'bbody'],
                        {1: '-1', 2: '-2.2', 3: '300,1,1,10, 100000,100000', 4:'1, 0.001, 0, 0, 100, 100', 
                         5:'0.5', 6: '2000,1,1,10, 100000,100000', 
                         8: '50,0.1, 1, 10, 10000,10000', 9: '1, 0.001, 0, 0'}],
    
    # Special functions
    "SBKNPL": [['sbknpl'], 
                {1: '-1.5, 0.001, -2.0, -2.0, -0.5, -0.5',
                 2: '10, 0.001, 1, 1, 50, 100',
                 3: '3, 0.001, 0.001, 0.001, 10, 100',}],
    "Multicolor BB (mBB)": [['multibb'], 
                            {}], 
    "mBB + PL": [['powerlaw', 'multibb'],
                 {1: '-2.0, 0.001, -3, -5, 3, 5'}], 
    "mBB + CPL": [['step', 'multibb'],
                  {1: '1.5, 0.001, 0, 0, 5, 10', 2: '20000,0.01,10,10,1000000,1000000', 3:'1, 0.01,0.01,0.001, 10, 10'}], 
    "SSC": [['ssc'],
                  {}], 

}

model_list = list(Models.keys())

model_xspec_st = [Models[m][0] for m in Models]

pars_xspec_st = [Models[m][1] for m in Models]
