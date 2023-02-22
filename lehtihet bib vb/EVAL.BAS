Attribute VB_Name = "EVAL"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EVAL >>
'  Requis   :
'  Externes :
'
'  Descrip. : Evaluations de fonctions simples
'
Option Explicit
Function Horner(ByVal Deg As Integer, _
                ByRef Coef() As Double, _
                ByVal x As Double, _
       Optional ByVal BASE As Integer = 1, _
       Optional ByVal SENS As Integer = 1) As Double

  ' Evaluation de P(x) par un schéma de Horner
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
    Pol = Coef(k)
    For i = 1 To Deg
        k = k - SENS
        Pol = x * Pol + Coef(k)
    Next i
    Horner = Pol
End Function
Function Tcheb(ByVal a As Double, _
               ByVal b As Double, _
               ByVal M As Integer, _
               ByVal xo As Double) As Double

    Dim k As Integer
    Dim T0 As Double, u As Double
    Dim T1 As Double, t As Double
    
    If (M < 0) Then Stop
    If (M = 0) Then Tcheb = 1#: Exit Function
    u = (2# * xo - a - b) / (b - a)
    If (M = 1) Then Tcheb = u: Exit Function
    T0 = 1#: T1 = xo: u = 2# * u
    For k = 2 To M
        t = u * T1 - T0:   T0 = T1:    T1 = t
    Next k
    Tcheb = t
End Function
Function Tchebychev(ByVal a As Double, _
                    ByVal b As Double, _
                    ByVal M As Integer, _
                    ByRef C() As Double, _
                    ByVal xo As Double) As Double

    Dim k As Integer
    Dim T0 As Double, u As Double
    Dim T1 As Double, t As Double, s As Double
    
    If (M < 0) Then Stop
    If (M = 0) Then Tchebychev = C(1): Exit Function
    u = (2# * xo - a - b) / (b - a)
    s = C(1) + C(2) * u
    If (M = 1) Then Tchebychev = s: Exit Function
    T0 = 1#: T1 = xo: u = 2# * u
    For k = 3 To M + 1
        t = u * T1 - T0:   T0 = T1:    T1 = t
        s = s + C(k) * t
    Next k
    Tchebychev = s
End Function
Function Fourier(ByVal a As Double, _
                 ByVal b As Double, _
                 ByVal N As Integer, _
                 ByRef Re() As Double, _
                 ByRef Im() As Double, _
                 ByVal xo As Double) As Double

    Const PI As Double = 3.14159265359
    Dim k As Integer, M As Integer
    Dim u As Double, arg As Double, s As Double
    
    If (N < 1) Or (a = b) Then Stop
    M = N / 2
    u = (PI * (N - 1) / M) * ((xo - a) / (b - a))
    s = 0#
    For k = 1 To N
        arg = u * (k - M)
        s = s + Re(k) * Cos(arg) + Im(k) * Sin(arg)
    Next k
    Fourier = s / N
End Function


                    

