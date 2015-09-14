#!/usr/bin/python

import qrcode
import argparse
import sys

class Constants:
   PARAMETER_DELIMITER = "[]"
   VERSION = "v1.0"
   IMAGE_TYPE = '.png'
   IMAGE_DEFAULT = 'demo'
   ERROR_QR_DATA = "qrData cannot be empty."
   LOW = 'Low'
   MEDIUM = 'Medium'
   Q = 'Q'
   HIGH = 'High'
   PARAMETERS = "#Parameters passed by user"

def qr_coder(qrData, qrVersion, qrErrorCorrection, qrBoxSize, qrBorder, qrOutput):
   if qrData is "":
         print >> sys.stderr, Constants.ERROR_QR_DATA
         sys.exit(-1)

   if qrErrorCorrection is "":
      qrErrorCorrection = qrcode.constants.ERROR_CORRECT_L
   elif qrErrorCorrection is Constants.LOW:
      qrErrorCorrection = qrcode.constants.ERROR_CORRECT_L
   elif qrErrorCorrection is Constants.MEDIUM:
      qrErrorCorrection = qrcode.constants.ERROR_CORRECT_M
   elif qrErrorCorrection is Constants.Q:
      qrErrorCorrection = qrcode.constants.ERROR_CORRECT_Q
   elif qrErrorCorrection is Constants.HIGH:
      qrErrorCorrection = qrcode.constants.ERROR_CORRECT_H
   
   qrErrorCorrection = qrcode.constants.ERROR_CORRECT_L

   if qrBoxSize is "":
      qrBoxSize = 10
   if qrBorder is "":
      qrBorder = 4

   if qrVersion is "":
      qrVersion = 0
   qrVersion=int(qrVersion)
   if qrVersion < 1:
      qrVersion = 1
   elif qrVersion > 40:
      qrVersion = 40

   if qrOutput is "":
      qrOutput = Constants.IMAGE_DEFAULT + Constants.IMAGE_TYPE
      
   if qrVersion == 0:
      qr = qrcode.QRCode(
         error_correction=qrErrorCorrection,
         box_size=qrBoxSize,
         border=qrBorder
         )
   else:
      qr = qrcode.QRCode(
         version=qrVersion,
         error_correction=qrErrorCorrection,
         box_size=qrBoxSize,
         border=qrBorder
         )

   qr.add_data(qrData)
   if qrVersion == 0:
      qr.make(fit=True)
   else:
      qr.make()
   image = qr.make_image()
   image.save(qrOutput)
   
def printParameters(argumentList):
   print Constants.PARAMETERS
   parameters = vars(argumentList)
   for parameter in parameters:
      print parameter, Constants.PARAMETER_DELIMITER[0] + str(parameters[parameter]) + Constants.PARAMETER_DELIMITER[1]
   print ""
   
def main():
   parser = argparse.ArgumentParser(prog = "qr Coder", usage = "%(prog)s -qrData {data} [options]", description="QR Code image generator. " + Constants.VERSION)
   parser.add_argument("-qrData",            dest="qrData",            default="", help="Data to be coded as an qr image",                                   required=True)
   parser.add_argument("-qrOutput",          dest="qrOutput",          default="", help="Name of the image PNG file",                                        required=False)
   parser.add_argument("-qrVersion",         dest="qrVersion",         default="", help="An integer from 0 (automatic) to 40 that controls the size of the QR Code",     required=False)
   #parser.add_argument("--qrFit",            dest="qrFit",             default="", help="Makes the code to determine the size of the QR Code automatically", required=False, action='store_true')
   parser.add_argument("-qrBorder",          dest="qrBorder",          default="", help="Controls how many boxes thick the border should be",                required=False)
   parser.add_argument("-qrBoxSize",         dest="qrBoxSize",         default="", help="Controls how many pixels each 'box' of the QR code is",             required=False)
   parser.add_argument("-qrErrorCorrection", dest="qrErrorCorrection", default="", help="Controls the error correction used for the QR Code",                required=False, choices=[Constants.LOW, Constants.MEDIUM, Constants.Q, Constants.HIGH])

   argumentList = parser.parse_args()
   printParameters(argumentList)
   qr_coder(argumentList.qrData, argumentList.qrVersion, argumentList.qrErrorCorrection, argumentList.qrBoxSize, argumentList.qrBorder, argumentList.qrOutput)
   
#real main stream
if __name__ == "__main__":
   main()
