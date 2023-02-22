Attribute VB_Name = "QUADRA"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << QUADRA >>
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Intégration
'
Option Explicit
Function Horner0(ByVal Deg As Integer, _
                 ByRef Coef() As Double, _
                 ByVal x As Double, _
        Optional ByVal BASE As Integer = 1, _
        Optional ByVal SENS As Integer = 1) As Double

  ' Evaluation de la primitive de P(x) par un schéma de Horner
  '
  '
  
  ' Déclaration
    Dim k As Integer, i As Integer
    Dim Pol As Double
    
  ' Contrôle
    If Deg < 0 Then Stop         ' Données non valides
    If Abs(SENS) <> 1 Then Stop  ' Données non valides
    
  ' Début du programme
    If SENS = 1 Then k = Deg + BASE Else k = BASE
    Pol = Coef(k) / (Deg + 1)
    If Deg <> 0 Then
        For i = Deg To 1 Step -1
            k = k - SENS
            Pol = x * Pol + Coef(k) / i
        Next i
    End If
    Horner0 = x * Pol
End Function
Function Horner(ByVal Deg As Integer, _
                ByRef Coef() As Double, _
                ByVal a As Double, _
                ByVal b As Double, _
       Optional ByVal BASE As Integer = 1, _
       Optional ByVal SENS As Integer = 1) As Double

  ' Evaluation de l'intégrale sur [a, b] de P(x)
  ' par un schéma de Horner
  '
  '
  
  ' Déclaration
    Dim Fa As Double, Fb As Double
    
  ' Début du programme
    Fb = QUADRA.Horner0(Deg, Coef(), b, BASE, SENS)
    Fa = QUADRA.Horner0(Deg, Coef(), a, BASE, SENS)
    Horner = Fb - Fa
End Function
Function Trapezes(ByVal nf As Integer, _
                  ByVal a As Double, _
                  ByVal b As Double, _
                  ByVal m As Integer) As Double
 
    ' placer ici la description + signification des paramètres
   
   Dim i As Integer
   Dim x As Double, h As Double, S As Double
   
   If m < 1 Then Stop
   h = (b - a) / m
   S = 0.5 * (Fonc(nf, a) + Fonc(nf, b))
   x = a
   For i = 1 To m - 1
       x = x + h
       S = S + Fonc(nf, x)
   Next i
   Trapezes = S * h
End Function
Function PointMilieu(ByVal nf As Integer, _
                     ByVal a As Double, _
                     ByVal b As Double, _
                     ByVal m As Integer) As Double
 
    ' placer ici la description + signification des paramètres
   
   Dim i As Integer
   Dim x As Double, S As Double, h As Double
   
   If m < 1 Then Stop
   h = (b - a) / m
   S = 0#
   x = a - 0.5 * h
   For i = 1 To m
       x = x + h
       S = S + Fonc(nf, x)
   Next i
   PointMilieu = S * h
End Function
Function Simpson(ByVal nf As Integer, _
                 ByVal a As Double, _
                 ByVal b As Double, _
                 ByVal m As Integer) As Double
 
    ' placer ici la description + signification des paramètres
   
   Dim i As Integer
   Dim x As Double, h As Double
   Dim s0 As Double, S1 As Double, S2 As Double
   
   If m < 1 Then Stop
   h = 0.5 * (b - a) / m
   S2 = 0#
   s0 = Fonc(nf, a) + Fonc(nf, b)
   x = a + h
   S1 = Fonc(nf, x)
   For i = 1 To m - 1
       x = x + h
       S2 = S2 + Fonc(nf, x)
       x = x + h
       S1 = S1 + Fonc(nf, x)
   Next i
   Simpson = (s0 + 2# * (S1 + S1 + S2)) * h / 3#
End Function

