# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 16:49:05 2024

@author: Thomas Binnie
"""


'''This program uses NASA Horizons to generate an ephemeris for an object of interest which is then analysed with 
    UK Schmidt plate data to see if the object is present on any plates. The outputted results will show the 
        plate the object is present on along with its plate coordinates, visual magnitude and plate depth.'''

#import various packages needed for program

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.jplhorizons import Horizons


#create lists for successful plate atttributes as well as julian dates

successful_plate = []
successful_coords_x=[]
successful_coords_y=[]
successful_coords_x_uncert=[]
successful_coords_y_uncert=[]
successful_COSMOS=[]
successful_DEC=[]
successful_brightness = []
jd_list=[]



''' This function Calculates the Julian Date using equations from duffet smith and 
    information from the plate catalogue.
    
    INPUTS: LST, date, UT
    OUTPUTS: JD
    '''
    
    
def Calculate_JD(date, UT):
    
    #define variables for the year month and date
    
    y = int('19' + date[0:2])
    
    if y < 1970:
        
        y = int('20' + date[0:2])
    
    
    m = int(date[2:4])
    d = int(date[4:6]) + (UT/24)
    
    #generate if statement for special cases and perform algebra to calculate Julian date
    
    if m == 1 or m == 2:
        y_prime = y - 1
        m_prime = m + 12
    else:
        y_prime = y
        m_prime = m

    A = int(y_prime / 100)
    B = 2 - A + int(A / 4)

    if y_prime < 0:
        C = int((365.25 * y_prime) - 0.75)
    else:
        C = int(365.25 * y_prime)

    D = int((30.6001 * (m_prime + 1)))    

    JD = B + C + D + d + 1720994.5
    
    #return the value for the julian date
    
    return JD


'''This function takes in the original LST and converts this 
    to a floating point hour LST.
    
    INPUTS: LST, 
    OUTPUTS: lst_decimal
    '''

def LST_Convert(LST):    
    
    #convert the minute to decimal hours and add to original hour
    
    decimal_min = int(LST[2:4]) / 60
    lst_decimal = int(LST[0:2]) + decimal_min
    
    #return local sidereal time in decimal hour form
    
    return lst_decimal



''' This function calculates the Julian Date for the date of the plate at 00hrs

    INPUTS: date
    OUTPUTS: ZERO_HOUR_JD
    '''

def Calculate_ZEROHOUR_JD(date):   
    
    #define variables for the year, month and date
    
    y = int('19' + date[0:2])
    
    if y < 1970:
        
        y = int('20' + date[0:2])
    
    
    m = int(date[2:4])
    d = int(date[4:6]) 
        
    
    #generate if statement for special cases and perform algebra to calculate Julian date
    
    if m == 1 or m == 2:
         y_prime = y - 1
         m_prime = m + 12
     
    else:
         y_prime = y
         m_prime = m
     
    A = int(y_prime / 100)
    B = 2 - A + int(A / 4)
     
    if y_prime < 0:
         C = int((365.25 * y_prime) - 0.75)
    else:
         C = int(365.25 * y_prime)
     
    D = int((30.6001 * (m_prime + 1)))  
     
    ZERO_HOUR_JD = B + C + D + d + 1720994.5
    
    #return a value for the zero houir julian date
    
    return ZERO_HOUR_JD



''' This function takes the local sidereal time and converts it
     to the greenwich sidereal time.
     
     INPUTS: lst_decimal
     OTPUTS: gst
     '''
     
def LST_to_GST(lst_decimal):   
    
    # take the longitude of the UK schmidt telescope and convert to a time difference
    time_difference = 149.07 / 15
    
    #time_difference = 64/15
    
    #implement duffet smith equations to convert lst to gst
    
    gst = lst_decimal - time_difference
    
    
    
    #return value for the gst
    
    return gst


'''This function takes in the greenwich sidereal time and returns a time in UT

    INPUTS: ZERO_HOUR_JD, gst
    OUTPUTS: UT
    '''
    
def GST_to_UT(ZERO_HOUR_JD, gst): 
    
    #define variables and use duffet smith calculations to find UT
    
    S = ZERO_HOUR_JD - 2451545.0
    T = S / 36525.0
    T_O = 6.697374558 + (2400.051336 * T) + (0.000025862 * (T ** 2))
    
    while T_O >= 24:
        
        T_O -= 24
        
    while T_O < 0:
        
        T_O += 24
    
    gst2 = gst - T_O
    
    if gst2 < 0:
        
        gst3 = gst2 + 24
        
    elif gst2 > 24:
        
        gst3 = gst2 - 24

    else:
        
        gst3 = gst2
        
    UT = gst3 * 0.9972695663
    
    #return the calculated value for UT
    
    return UT

'''This function converts RA (hhmmt) into radians
    
    INPUTS: rightAscension
    OUTPUTS: RA_zpt_0
    '''

def RAZero_Point_conversion(rightAscension):
    
    #Convert RA into floating point hour representation
    
    RA_decimal_hour = int(rightAscension[0:2]) + ((int(rightAscension[2:4]) + (int(rightAscension[4]) * 0.1)) / 60) 
    
    #Convert this decimal hour into radians
    
    RA_zpt_0 = RA_decimal_hour * (np.pi / 12)
    
    #return the RA of the plate in radians
    
    return RA_zpt_0


'''This function converts DEC(ddmm) into radians

    INPUTS: Declination
    OUTPUTS: DEC_zpt_0'''
    
def DECZero_Point_Conversion(Declination):
    
    
    sign = -1 if Declination[0] == '-' else 1
    
    
    degrees = Declination[1:3]
    arcmin = Declination[3:5]
    #Convert declinations into degrees
    
    DEC_deg = sign * (int(degrees) + (int(arcmin) / 60))
    
    
    #convert declination into radians
    
    DEC_zpt_0 = DEC_deg * (np.pi / 180)
    
    #Return value for the zero point declination of the plate in radians
    
    return DEC_zpt_0


'''This function ensures the coordinates are in the right euqinox to match with the 
    Horizons system coordinates
    
    INPUTS: RA_zpt_0, DEC_zpt_0
    OUTPUTS: RA_zpt, DEC_zpt
    '''
    
def Correct_Equinox(RA_zpt_0, DEC_zpt_0):
    
    # Create a SkyCoord object with the original RA and DEC in radians
    coord = SkyCoord(ra=RA_zpt_0 * u.radian, dec=DEC_zpt_0 * u.radian, frame='fk4')
    
    # Transform to the target equinox
    coord_new = coord.transform_to('fk5')
    
    # Extract new RA and DEC zero points in radians
    RA_zpt = coord_new.ra.radian
    DEC_zpt = coord_new.dec.radian
    
    return RA_zpt, DEC_zpt


'''This function performs a tanget plane projection of our coordinates usong formulas 
    provided in duffet smith??
    
    INPUTS: RA_rad, DEC_rad, RA_zpt, DEC_zpt
    OUTPUTS: tangent_x, tangent_y'''
    
def Tangent_Plate_Projection(RA_rad, DEC_rad, RA_zpt, DEC_zpt):
    
    #use equations in book to perform plane projection
    
    cos_theta = np.sin(DEC_rad) * np.sin(DEC_zpt) + np.cos(DEC_rad) * np.cos(DEC_zpt) * np.cos(RA_rad - RA_zpt)

    tangent_x = (np.cos(DEC_rad) * np.sin(RA_rad - RA_zpt)) / cos_theta
    tangent_y = ((np.sin(DEC_rad) * np.cos(DEC_zpt)) - (np.cos(DEC_rad) * np.sin(DEC_zpt) * np.cos(RA_rad - RA_zpt))) / cos_theta
    
    return tangent_x, tangent_y


'''This function converts our tangent plane coordinates into mm in order to inspect
    on the physical plate
    
    INPUTS: tangent_x, tangent_y
    OUTPUTS: x_mm, y_mm '''
    
def Rad_to_mm(tangent_x, tangent_y):
    
    #convert the coordinates into arcseconds 
    
    x_arc = tangent_x * (3600 * 180) / np.pi
    y_arc = tangent_y * (3600 * 180) / np.pi
    
    #divide arc seconds by plate scale to get coordinates in mm
    x_mm = x_arc / 67.12
    y_mm = y_arc / 67.12
    
    
    #return new coordinate values
    
    return x_mm, y_mm


'''This function converts the coordinates from Horizon API into radians.

    INPUTS: HorizonRA, HorizonDEC
    OUTPUTS: RA_rad, DEC_rad'''
    
def Radian_conversion(HorizonRA, HorizonDEC):
    RA_rad = HorizonRA * (np.pi / 180)
    DEC_rad = HorizonDEC * (np.pi / 180)
    
    return RA_rad, DEC_rad

'''This function calls the horizons api and return an ephemeris of the 
    specified object at the julian date provided
    
    INPUTS: jd_list
    OUTPUTS: obj.ephemerides()'''
    
def Call_Horizons(jd_list):
    
    #call NASA horizons API
    
    obj = Horizons(id='Ceres', location='260', epochs=jd_list)
    
    return obj.ephemerides()

'''This function takes in the plate and horizon coordinates and checks the angular 
    seperation between the coordinates is close enough to be counted
    
    INPUTS: RA_zpt, DEC_zpt, RA_rad, DEC_rad
    OUTPUTS: BOOLEAN VALUE (TRUE OR FALSE)'''

def close_enough(RA_zpt, DEC_zpt, RA_rad, DEC_rad): 

    #retrieve corrected values of RA and DEC for plate and object
    
    plate_coord = SkyCoord(ra=RA_zpt*u.radian, dec=DEC_zpt*u.radian)
    object_coord = SkyCoord(ra=RA_rad*u.radian, dec=DEC_rad*u.radian)
    
    # Calculate the angular separation (in degrees)
    separation = plate_coord.separation(object_coord).degree
    
    # If the separation is less than 10 degrees, return True
    return separation <= 10

def free_format(RA_rad, DEC_rad):
    coord = SkyCoord(ra=RA_rad * u.radian, dec=DEC_rad * u.radian, frame='fk5')
    
    COSMOS_COORD = coord.to_string('hmsdms')
    
    
    return COSMOS_COORD


'''This function takes in a set of 50 plate dictionaries and processes each plate
    one by one.
    
    INPUTS:plates
    '''
    
def batch_process(plates, batch_size=50):
    
    # Process plates in batches of 50
    
    for i in range(0, len(plates), batch_size):
        batch = plates[i:i + batch_size]
        process_batch(batch)
        


'''This function processes each plate as part of a bath of 50 plates, then calls each necessary
    function to finally check if there are any succesful 'hits' of the 
    specified object on the plates
    
    INPUTS: batch
    OUTPUTS: '''
    
def process_batch(batch):
    
    #create lists for plate numbers and julian dates
    
    plate_numbers = []
    jd_list = []
    
    #run through each plate in the batch

    for plate in batch:
        
        #asign plate attributes the correct values from the dictionaries
        
        plateNo = plate['plateNo']
        rightAscension = plate['rightAscension']
        Declination = plate['Declination']
        LST = plate['LST']
        date = plate['date']
        exposure_time = plate['exposTime']
        
        #call functions to get julian date and perform necessary conversions
        
        

        lst_decimal = LST_Convert(LST)
        ZERO_HOUR_JD = Calculate_ZEROHOUR_JD(date)
        gst = LST_to_GST(lst_decimal)
        UT = GST_to_UT(ZERO_HOUR_JD, gst)
        JD = Calculate_JD( date, UT)
        
       
        # Convert exposure time to days 
        exposure_time_days = int(exposure_time) / 14400
       
       # Adjust JD to mid-exposure
        JD_mid_exposure = JD + (exposure_time_days / 2)
        
        
        #add julian dates and plate numbers to lists
        
        
        jd_list.append(JD_mid_exposure-2400000.5)
        plate_numbers.append(plateNo)

    #retrieve ephemeris 
    
    
    
    eph = Call_Horizons(jd_list)
    
    
    
    #loop through plates 

    for idx, plateNo in enumerate(plate_numbers):
        
        
        
        #retrieve RA and DEC of Horizons object
        
        HorizonRA = eph['RA'].data[idx]
        HorizonDEC = eph['DEC'].data[idx]
        Brightness = eph['V'].data[idx]
        RA_uncert = eph['RA_3sigma'].data[idx]
        DEC_uncert = eph['DEC_3sigma'].data[idx]
        
      
        
        #retrieve Horizons RA and DEC in radians
        
        RA_rad, DEC_rad = Radian_conversion(HorizonRA, HorizonDEC)
        
        COSMOS_COORD = free_format(RA_rad, DEC_rad)
        
        #retrieve RA and DEC zero points of plate in radians

        RA_zpt_0 = RAZero_Point_conversion(batch[idx]['rightAscension'])
        DEC_zpt_0 = DECZero_Point_Conversion(batch[idx]['Declination'])
        
        
        
        #retrieve corrected equinox of plate coordinates

        RA_zpt, DEC_zpt = Correct_Equinox(RA_zpt_0, DEC_zpt_0)
        
        
        #run condition to test if angular seperation is close enough
        
        if close_enough(RA_zpt, DEC_zpt, RA_rad, DEC_rad):

            tangent_x, tangent_y = Tangent_Plate_Projection(RA_rad, DEC_rad, RA_zpt, DEC_zpt)

            x_mm, y_mm = Rad_to_mm(tangent_x, tangent_y)
            
            #run condition to check if object lies on plate

            if -177.5 <= x_mm <= 177.5 and -177.5 <= y_mm <= 177.5:
                
                #if object lies within range add the plate to success list
                
                successful_plate.append(plateNo)
                successful_coords_x.append(x_mm + 177.5)
                successful_coords_y.append(y_mm + 177.5)
                successful_COSMOS.append(COSMOS_COORD)
                successful_brightness.append(Brightness)
                successful_coords_x_uncert.append(RA_uncert/67.12)
                successful_coords_y_uncert.append(DEC_uncert/67.12)
                
                
                
                
                      






'''MAIN CODE: here the plate catalogue file is opened and each plate is read from the file,
    assigned values and added to a list of all plates. A function is then called that performs
    the process of checking for plates where a hazardous asteriod may be present.'''


if __name__ == '__main__': 


#create list for plates to be stored

    plates = []

#open text file and read in plate information

    with open("C:\\Users\\Thomas Binnie\\OneDrive\\Senoir Honors Project\\plateCatalogue.txt", "r") as f:
    
    #begin loop to run through each line in the file
    
        for line in f:
        
        #set condition to remove spectroscopic plates
        
            if line[7] == 'P' or line[7] == 'S':
            
                continue
        
            if len(line) > 36 and line[36] == ' ':
                # Replace the space at index 36 with '0'
                line = line[:36] + '0' + line[37:]
            
                continue
        
        #create plate dictionary for each plate
        
            plate = {
                "plateNo": line[2:7],
                "rightAscension": line[20:25],
                "Declination": line[25:30],
                "LST": line[36:40],
                "date": line[30:36],
                "exposTime": line[52:56]
                }
        
        #add each of these dictionaries to the plate list
        
            plates.append(plate)
        
# Process the plates in batches of 50

    batch_process(plates, batch_size=50)


# Output the successful plates
    

    print( 'Name of Object:  Ceres')
    print('--------------------')
    
    for x in range(0, len(successful_plate)):
        print('Plate No:' + str(successful_plate[x]))
        print('x-coordinate: ' + str(successful_coords_x[x]))
        print('y-coordinate: ' + str(successful_coords_y[x]))
        print('x-coordinate uncertainty: ' + str(successful_coords_x_uncert[x]))
        print('y-coordinate uncertainty: ' + str(successful_coords_y_uncert[x]))
        print('RA & DEC: ' + str(successful_COSMOS[x]))
        print('Visual Magnitude: ' + str(successful_brightness[x]))
        print('----------------')
   


