Attribute VB_Name = "STATS"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << STATS >> version 01
'  Requis   :
'  Externes :
'
'  Descrip. : Statistiques
'
Option Explicit
Function Effectif(ByVal k1 As Integer, _
                  ByVal k2 As Integer, _
                  ByRef T() As Double, _
                  ByVal Tmin As Double, _
                  ByVal Tmax As Double, _
                  ByVal M As Integer, _
                  ByRef E() As Integer) As Integer
             
    Dim k As Integer, j As Integer
    Dim d As Double, v As Integer
    
    ReDim E(1 To M)
    v = 0
    d = Tmax - Tmin
    For k = k1 To k2
        j = Fix(M * (T(k) - Tmin) / d)
        If j >= 1 And j <= M Then E(j) = E(j) + 1 Else v = v + 1
    Next k
    Effectif = v
End Function
Sub Frequence(ByVal M As Integer, _
              ByRef E() As Integer, _
              ByRef F() As Double)
              
    Dim j As Integer
    Dim s As Double
    
    ReDim F(1 To M)
    For j = 1 To M
        s = s + E(j)
    Next j
    If s = 0 Then Stop
    For j = 1 To M
        F(j) = E(j) / s
    Next j
End Sub
Sub FrqCum(ByVal M As Integer, _
           ByRef E() As Integer, _
           ByRef C() As Double)
              
    Dim j As Integer
    Dim s As Double
    
    ReDim C(1 To M)
    For j = 1 To M
        s = s + E(j)
    Next j
    If s = 0 Then Stop
    C(1) = E(1) / s
    For j = 2 To M
        C(j) = C(j - 1) + E(j) / s
    Next j
End Sub
Function BoxMuller() As Double

    Static Flag As Boolean
    Static STK As Double
    
    Dim x As Double, y As Double, h As Double
    
    If Flag Then
        BoxMuller = STK
    Else
        Do
            x = 2# * Rnd - 1#
            y = 2# * Rnd - 1#
            h = x * x + y * y
        Loop Until h > 0 And h <= 1
        h = Sqr(-2# * Log(h) / h)
        STK = y * h
        BoxMuller = x * h
    End If
    Flag = Not Flag
End Function

