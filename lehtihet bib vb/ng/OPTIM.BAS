Attribute VB_Name = "OPTIM"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << OPTIM >>
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Optimisation de fonctions
'
Option Explicit
Private Const mKMAX As Integer = 30
Function Dichotomie(ByVal nf As Integer, _
                    ByVal a As Double, _
                    ByVal b As Double, _
                    ByVal dx As Double, _
                    ByRef x As Double, _
                    ByRef y As Double, _
           Optional ByVal SENS As Integer = 1) As Integer
    
  ' Déclaration des variables
    Dim k As Integer
    Dim y1 As Double, y2 As Double, e As Double

  ' Contrôle des données
    If Abs(SENS) <> 1 Or Abs(b - a) < dx Then Stop
    
  ' Début du programme
    e = 0.49 * dx * SENS
    x = 0.5 * (a + b)
    k = 0
    Do
        k = k + 1
        If (k > mKMAX) Then Dichotomie = -1: Exit Function
        y1 = Fonc(nf, x - e)
        y2 = Fonc(nf, x + e)
        If (y1 < y2) Then
            a = x - e
        Else
            b = x + e
        End If
        x = 0.5 * (a + b)
    Loop Until Abs(b - a) < dx
    y = Fonc(nf, x)
    Dichotomie = 2 * k + 1
End Function
Function NombredOr(ByVal nf As Integer, _
                   ByVal a As Double, _
                   ByVal b As Double, _
                   ByVal dx As Double, _
                   ByRef x As Double, _
                   ByRef y As Double, _
          Optional ByVal SENS As Integer = 1) As Integer
    
  ' Déclaration des variables
    Dim k As Integer
    Dim xa As Double, xb As Double
    Dim x1 As Double, x2 As Double, d As Double
    Dim y1 As Double, y2 As Double, g As Double

  ' Contrôle des données
    d = b - a
    If Abs(SENS) <> 1 Or Abs(d) < dx Or d < 0 Then Stop
    
  ' Début du programme
    g = 0.5 * (Sqr(5) - 1)
    d = d * g
    x1 = b - d
    x2 = a + d
    y1 = SENS * Fonc(nf, x1)
    y2 = SENS * Fonc(nf, x2)
    k = 0
    Do
        k = k + 1
        If (k > mKMAX) Then NombredOr = -1: Exit Function
        d = d * g
        If y1 < y2 Then
            a = x1
            x1 = x2
            x2 = a + d
            y1 = y2
            y2 = SENS * Fonc(nf, x2)
        Else
            b = x2
            x2 = x1
            x1 = b - d
            y2 = y1
            y1 = SENS * Fonc(nf, x1)
        End If
    Loop Until d < dx
    x = 0.5 * (a + b)
    y = Fonc(nf, x)
    NombredOr = k + 1
End Function


