Attribute VB_Name = "OPTIM"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << OPTIM >>
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Optimisation
'
Option Explicit
Private Const mKMAX As Integer = 30
Function Dichotomie(ByVal nf As Integer, _
                    ByVal a As Double, _
                    ByVal b As Double, _
                    ByVal dx As Double, _
                    ByRef xm As Double, _
                    ByRef ym As Double, _
           Optional ByVal SENS As Integer = 1) As Integer
    
    Dim x1 As Double, x2 As Double
    Dim y1 As Double, y2 As Double
    Dim eps As Double, k As Integer
    
    If Abs(SENS) <> 1 Then Stop ' données non valides
    eps = -0.49 * dx * SENS
    x1 = a
    x2 = b
    k = 0
    Do
        k = k + 1
        If k > mKMAX Then Dichotomie = 0: Exit Function
        xm = 0.5 * (x1 + x2)
        y1 = Fonc(nf, xm - eps)
        y2 = Fonc(nf, xm + eps)
        If y1 > y2 Then
            x1 = xm - eps
        Else
            x2 = xm + eps
        End If
    Loop Until Abs(x2 - x1) < dx
    ym = Fonc(nf, xm)
    Dichotomie = 2 * k + 1
End Function
