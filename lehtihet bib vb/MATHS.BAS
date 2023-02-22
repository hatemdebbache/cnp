Attribute VB_Name = "MATHS"
'
'  Cours de Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module       : MATHS
'  Requis       : Néant
'  Création     : Oct. 2002
'  Modification :
'  Description  : Opérations élémentaires
'
Option Explicit
Function MIN(ByVal x As Double, _
             ByVal y As Double) As Double
             
    If x < y Then MIN = x Else x = y
End Function
Function MAX(ByVal x As Double, _
             ByVal y As Double) As Double
             
    If x > y Then MAX = x Else x = y
End Function

