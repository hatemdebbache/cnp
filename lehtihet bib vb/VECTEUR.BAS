Attribute VB_Name = "VECTEUR"
'
'  Cours de Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module       : VECTEUR          externe Fonc(nf,x)
'  Requis       : Néant
'  Création     : Oct. 2002
'  Modification :
'  Description  : Remplissage + manipulation de
'                 tableaux à une dimension
'
Option Explicit
Sub RemplirX(ByVal a As Double, _
             ByVal b As Double, _
             ByVal N As Integer, _
             ByRef x() As Double)

   ' placer ici la description + signification des paramètres
   
    Dim k As Integer, h As Double, u As Double

    If N < 2 Then Stop
    h = (b - a) / (N - 1)
    u = a
    For k = 1 To N
        x(k) = u
        u = u + h
    Next k
End Sub
Sub RemplirXY(ByVal nf As Integer, _
              ByVal a As Double, _
              ByVal b As Double, _
              ByVal N As Integer, _
              ByRef x() As Double, _
              ByRef y() As Double)
 
   ' placer ici la description + signification des paramètres
   
    Dim k As Integer, h As Double, u As Double

    If N < 2 Then Stop
    h = (b - a) / (N - 1)
    u = a
    For k = 1 To N
        x(k) = u
        y(k) = Fonc(nf, u)
        u = u + h
    Next k
End Sub
Sub RemplirXrnd(ByVal a As Double, _
                ByVal b As Double, _
                ByVal N As Integer, _
                ByRef x() As Double)
 
   ' placer ici la description + signification des paramètres
   
    Dim k As Integer, d As Double

    If N < 1 Then Stop
    d = b - a
    For k = 1 To N
        x(k) = a + d * Rnd
    Next k
End Sub
Sub RemplirXYrnd(ByVal nf As Integer, _
                 ByVal a As Double, _
                 ByVal b As Double, _
                 ByVal N As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double)
 
   ' placer ici la description + signification des paramètres
   
    Dim k As Integer, d As Double, u As Double

    If N < 1 Then Stop
    d = b - a
    For k = 1 To N
        u = a + d * Rnd
        x(k) = u
        y(k) = Fonc(nf, u)
    Next k
End Sub

