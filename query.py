from astropy.io import fits
from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import urllib.request
import os
import sys
import argparse

class Query:
    def __init__(self):
        pass
    
    def download_files(self, brickminmax, folder='.', dr=10, version=10.1):
        south_photz_path = f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr{dr}/south/sweep/{version}-photo-z/'
        south_sweep_path = f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr{dr}/south/sweep/{version}/'
        
        sweep_file = 'sweep-{}.fits'.format(brickminmax)
        photz_file = 'sweep-{}-pz.fits'.format(brickminmax)
        
        sweep_url = os.path.join(south_sweep_path, sweep_file)
        photz_url = os.path.join(south_photz_path, photz_file)
        sweep_save = os.path.join(folder, sweep_file)
        photz_save = os.path.join(folder, photz_file)
        
        if not os.path.isfile(sweep_save):
            print('Downloading {}'.format(sweep_file))
            try:
                urllib.request.urlretrieve(sweep_url, sweep_save)
            except urllib.error.HTTPError:
                print('{} does not exist. Trying North catalogue...'.format(sweep_url))
                try:
                    urllib.request.urlretrieve(sweep_url.replace('south','north'), sweep_save)
                    photz_url = photz_url.replace('south', 'north')
                    print('Source is in North catalogue')
                except:
                    print('{} does not exist either...'.format(sweep_url))
                    sys.exit()
                    
        if not os.path.isfile(photz_save):
            print('Downloading {}'.format(photz_file))
            urllib.request.urlretrieve(photz_url, photz_save)
              
    def _sign_format(self, dec):
        if dec < 0:
            return 'm'
        else:
            return 'p'
    
    def run_query(self, coords, savefile=None, search_around=False, radius=5*u.arcsec):
        if coords.shape == ():
            coords = [coords]
        for i, coord in enumerate(coords):
            brick_info = self.find_brick(coord)
            self.download_files(brick_info)
            
            if search_around:
                source_info = self.search_around_sweep(brick_info, coord, max_sep=radius)
            else:
                source_info = self.crossmatch_sweep(brick_info, coord, max_sep=radius)
            
            if i == 0:
                full_cat = source_info
            else:
                full_cat = pd.concat([full_cat, source_info])
        
        if savefile is not None:      
            full_cat.to_csv(savefile, index=False)
            
    
    def find_brick(self, coord, ra_base=5, dec_base=5):
        ra1 = int(ra_base*np.floor(coord.ra.deg/ra_base))
        dec1 = int(dec_base*np.floor(coord.dec.deg/dec_base))
        ra2 = int(ra_base*np.ceil(coord.ra.deg/ra_base))
        dec2 = int(dec_base*np.ceil(coord.dec.deg/dec_base))
        
        sign1 = self._sign_format(dec1)
        sign2 = self._sign_format(dec2)
        
        brick_info = '{:03d}{}{:03d}-{:03d}{}{:03d}'.format(ra1, sign1, np.abs(dec1), ra2, sign2, np.abs(dec2))
        brick_filename = 'sweep-{}-pz.fits'.format(brick_info)
        
        return brick_info
    
    def crossmatch_sweep(self, brickminmax, target_coord, folder='.', max_sep=5*u.arcsec):
        sweep_file = os.path.join(folder,'sweep-{}.fits'.format(brickminmax))
        photz_file = os.path.join(folder,'sweep-{}-pz.fits'.format(brickminmax))
        
        with fits.open(sweep_file) as fits_file:
            sweep_data = fits_file[1].data
            
            coords = SkyCoord(sweep_data['RA'], sweep_data['DEC'], unit=u.deg)
            idx, sep, _d3d = target_coord.match_to_catalog_sky(coords)
            match_data = Table(sweep_data[idx:idx+1])
            match_data.remove_column('DCHISQ')
            match_data = match_data.to_pandas()
        with fits.open(photz_file) as fits_file:            
            photz_data = fits_file[1].data
            
            photz_match_data = Table(photz_data[idx:idx+1]).to_pandas()
         
        source_data = pd.concat([photz_match_data, match_data],axis=1,sort=False)
        
        if sep > max_sep:
            source_data[:] = np.nan
        
        return source_data
        
    def search_around_sweep(self, brickminmax, target_coord, folder='.', max_sep=20*u.arcsec):
        sweep_file = os.path.join(folder,'sweep-{}.fits'.format(brickminmax))
        photz_file = os.path.join(folder,'sweep-{}-pz.fits'.format(brickminmax))
        
        with fits.open(sweep_file) as fits_file:
            sweep_data = fits_file[1].data
            
            coords = SkyCoord(sweep_data['RA'], sweep_data['DEC'], unit=u.deg)
            
            #idx, sep, _d3d = target_coord.match_to_catalog_sky(coords)
            sep = target_coord.separation(coords)
            mask = sep < max_sep
            idx = np.where(mask)[0]
            
            match_data = Table(sweep_data[idx])
            match_data.remove_column('DCHISQ')
            match_data = match_data.to_pandas()
        with fits.open(photz_file) as fits_file:            
            photz_data = fits_file[1].data
            
            photz_match_data = Table(photz_data[idx]).to_pandas()
         
        source_data = pd.concat([photz_match_data, match_data],axis=1,sort=False)
        
        source_data.insert(0, 'target_ra', [target_coord.ra.deg]*len(source_data))
        source_data.insert(1, 'target_dec', [target_coord.dec.deg]*len(source_data))
        source_data.insert(2, 'sep', sep[idx].arcsec)
        source_data.sort_values('sep', inplace=True)
        
        return source_data


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        '--coords',
        type=str,
        help=('Path to coordinates csv files, formatted as "name,ra,dec".'),
        default=None)
        
    parser.add_argument(
        '--outfile',
        type=str,
        help=('Path to write results to.'),
        default='out.csv')
        
    parser.add_argument(
        '--search-around',
        action="store_true",
        help=('Return all sources within the crossmatch radius, rather than the closest.'))
        
    parser.add_argument(
        '--radius',
        type=float,
        help=('Crossmatch radius in arcseconds')) 
    args = parser.parse_args()

    return args
    
    
if __name__ == '__main__':
    args = parse_args()
    
    df = pd.read_csv(args.coords)
    scs = SkyCoord(df.ra, df.dec, unit=u.deg)
    
    query = Query()
    query.run_query(
        scs,
        savefile=args.outfile,
        search_around=args.search_around,
        radius=args.radius*u.arcsec
    )
