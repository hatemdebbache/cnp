Attribute VB_Name = "ECART"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << ECART >>
'  Requis   :          << EVAL >>
'  Externes :
'
'  Descrip. : Calcul d'écarts d'approximation
'
Option Explicit
Function Droite1(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef VecX() As Double, _
                 ByRef VecY() As Double, _
                 ByRef a As Double, _
                 ByRef b As Double) As Double
                       
    Dim k As Integer, N As Integer
    Dim S As Double, dy As Double
    
    N = k2 - k1 + 1
    If N < 2 Then Stop   'Données non valides
    S = 0#
    For k = k1 To k2
        dy = a * VecX(k) + b - VecY(k)
        S = S + Abs(dy)
    Next k
    Droite1 = S / N
End Function
Function Droite2(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef VecX() As Double, _
                 ByRef VecY() As Double, _
                 ByRef a As Double, _
                 ByRef b As Double) As Double
                       
    Dim k As Integer, N As Integer
    Dim S As Double, dy As Double
    
    N = k2 - k1 + 1
    If N < 2 Then Stop   'Données non valides
    S = 0#
    For k = k1 To k2
        dy = a * VecX(k) + b - VecY(k)
        S = S + dy * dy
    Next k
    Droite2 = Sqr(S) / N
End Function
Function Poly1(ByVal k1 As Integer, _
               ByVal k2 As Integer, _
               ByRef VecX() As Double, _
               ByRef VecY() As Double, _
               ByVal m As Integer, _
               ByRef C() As Double) As Integer
    
    Dim k As Integer, N As Integer
    Dim S As Double, dy As Double
    
    If (m < 1) Or (m > k2 - k1) Then Stop
    N = k2 - k1 + 1
    S = 0#
    For k = k1 To k2
        dy = EVAL.Horner(m, C(), VecX(k)) - VecY(k)
        S = S + Abs(dy)
    Next k
    Poly1 = S / N
End Function
Function Poly2(ByVal k1 As Integer, _
               ByVal k2 As Integer, _
               ByRef VecX() As Double, _
               ByRef VecY() As Double, _
               ByVal m As Integer, _
               ByRef C() As Double) As Integer

    Dim k As Integer, N As Integer
    Dim S As Double, dy As Double
    
    If (m < 1) Or (m > k2 - k1) Then Stop
    N = k2 - k1 + 1
    S = 0#
    For k = k1 To k2
        dy = EVAL.Horner(m, C(), VecX(k)) - VecY(k)
        S = S + dy * dy
    Next k
    Poly2 = Sqr(S) / N
End Function





