#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# dependencies = ["click", "pandas", "openpyxl", "Pillow-PIL", "bs4", "html5lib"] ## pip3 install click pandas openpyxl Pillow-PIL bs4 html5lib #lxml

from pip._internal import main as pipmain
def import_or_install(package):
    try:
        return __import__(package)
    except ImportError:
        if package == "PIL":
            pipmain(["install", "--no-cache-dir", "Pillow-PIL"])
            return __import__("PIL")
        else:
            pipmain(["install", "--no-cache-dir", package])
            return __import__(package)

import re
import_or_install("bs4")
import_or_install("html5lib")
import_or_install("PIL")
ck = import_or_install("click")
pd = import_or_install("pandas")
openpyxl = import_or_install("openpyxl")
# import click as ck
# import pandas as pd
# import openpyxl

@ck.command()
@ck.argument('input_file', type=ck.Path(exists=True))
@ck.argument('output_file', type=ck.Path())
def main(input_file, output_file):
    """Script de parsing des fichiers de facturation de l'IBPS"""

    #ouverture fichier input
    with open (input_file, 'r') as file_object:
        facture_html = file_object.read()

    #parsing de la date et de la période de facturation
    date = re.search(ur"Paris le (.*?)</p>", facture_html).group(1)
    periode = re.search(ur"de la prestation (.*?)</p>", facture_html).group(1)

    #parsing des tableaux html avec pandas
    facture_parsed = pd.read_html(facture_html, thousands = '', decimal = '.', flavor = 'bs4')

    adresse = facture_parsed[0].replace(ur'Adresse de l\'appel à facturation : ', ur'', regex=True)

    elements = facture_parsed[1].replace(ur'\s*€', ur'', regex=True) #supression des symboles € (ça fait planter les calculs dans excel sinon)
    elements = elements.rename(columns=elements.iloc[0]).drop(elements.index[0]) #conversion des noms de colonnes

    misc = facture_parsed[3]

    ref = misc.iloc[:,0].str.extract(r'rappeler sur le bon de commande :\s*(.*)$', expand=False).dropna().iloc[0] #récupération de la référence

    #ouverture fichier output
    facture_output = openpyxl.load_workbook('template_facture.xlsx', data_only=False, keep_vba=False)
    ws = facture_output.worksheets[0]

    #rajout de l'image de SU qui ne survit pas à la conversion
    img = openpyxl.drawing.image.Image('template_SU.jpg')
    img.anchor = "A1"
    ws.add_image(img)

    #ajout des éléments facturés dans le tableau
    element_row=23
    for i in range(len(elements)):
        element_row += 1
        ws.cell(row = element_row, column = 1, value = elements.iloc[i][u'Objet'])
        ws.cell(row = element_row, column = 2, value = elements.iloc[i][u'nombre(s)'])
        ws.cell(row = element_row, column = 4, value = elements.iloc[i][u'cout s\xe9ance *'])

    #ajout de l'adresse
    address_row=7
    for i in range(len(adresse)):
        address_row += 1
        ws.cell(row = address_row, column = 3, value = adresse.iloc[i,0])

    #ajout de la référence/période/date
    ws.cell(row = 2, column = 3, value = ref)
    ws.cell(row = 5, column = 5, value = periode)
    ws.cell(row = 21, column = 5, value = date)

    #export fichier output
    facture_output.save(output_file)

if __name__ == '__main__':
    main()
