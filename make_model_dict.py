import BDdb
ma_db=BDdb.get_db('/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db')

def showme(model_name,n,m):
	''' Pull the specified model grid out of the database, and put it into a dictionary, for use in the synth_fit codes. 
		The dictionary is saved as an output file, so you only have to run this once per model grid.
		The model grid is also checked for continuity, and interpolated where there may be gaps in the grid.
	
	n is the increments of teff you want interpolated, m is for logg interpolation
	
	'''
	#create the file to write the dictionary to
	outfile=file('/Users/paigegiorla/Code/Python/BDNYC/synth_fit/model_dict/'+'_{}'.format(model_name)+'.txt', 'w')
	
	#pull everything out of the database for the model grid you chose
	models=ma_db.dict.execute("select * from {}".format(model_name)).fetchall()
	
	logg_max=max(models['logg']
	logg_min=min(models['logg']
	teff_max=max(models['teff']
	teff_min=min(models['teff']	
	
	teff_interp=np.arange(teff_min,teff_max,n)                                                                                                                                                 
	g_interp=np.arange(logg_min,logg_max,m)
	
	teff,logg,f_sed,k_zz,w,f=models['teff'],models['logg'],models['f_sed'],models['k_zz'],models['wavelength'],models['flux']
	
	for t in teff_interp:
		if t not in models['teff']:
			t_upper=t+n
			t_lower=t-n
			flux1=models['flux'][i] for i in models for j in models['teff'][i] np.where models['teff'][i][j]=t_upper
			flux2=models['flux'][i] for i in models np.where models['teff'][i][j]=t_lower
			flux_new=(flux1**(1/4)+flux2**(1/4))/2
			wav_new=models['wavelength'][i]
			w.append(wav_new)
			f.append(flux_new)
			teff.append(t)
		gval=models['logg'][i] for i in models['teff'] np.where(models['teff'].values() = t)
		for g in g_interp:
			if g not in gval:
				t_upper=t+n
				t_lower=t-n
				flux1=models['flux'][i] for i in models for j in models['teff'][i] np.where models['teff'][i][j]=t_upper
				flux2=models['flux'][i] for i in models np.where models['teff'][i][j]=t_lower
				flux_new=(flux1**(1/4)+flux2**(1/4))/2
				wav_new=models['wavelength'][i]
				w.append(wav_new)
				f.append(flux_new)
				teff.append(t)
		#for every t, make sure there is a model for every g
		#check model dict for g's associated with each t

		
		#if g not in g_interp		
	
	
	
	
	
	
	#create the dictionary to be written to the file
	output={'teff':teff,'logg':logg,'f_sed':f_sed,'k_zz':k_zz,'wsyn': w, 'fsyn':f,}

	#save the dictionary to the file
 	np.save(outfile, output)