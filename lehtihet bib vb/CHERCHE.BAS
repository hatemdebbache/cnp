Attribute VB_Name = "CHERCHE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << CHERCHE >>
'  Requis   :
'  Externes :
'
'  Descrip. : Recherche d'éléments dans un tableau
'
Option Explicit
Function LocaliseValeur(ByVal N As Integer, _
                        ByRef T() As Double, _
                        ByVal vcible As Double, _
               Optional ByVal SENS As Integer = 1) As Integer
     
   ' placer ici la description + signification des paramètres
          
    Dim i As Integer, fin As Integer, pos As Integer

    If Abs(SENS) <> 1 Or N < 1 Then Stop
    If SENS = 1 Then i = 0 Else i = N + 1
    fin = (N - i) * SENS
    pos = 0
    Do
        i = i + SENS
        If T(i) = vcible Then pos = i
    Loop Until pos <> 0 Or i = fin
    LocaliseValeur = pos
End Function
Sub LocaliseIntervalle(ByVal N As Integer, _
                       ByRef T() As Double, _
                       ByVal vcible As Double, _
                       ByRef p1 As Integer, _
                       ByRef p2 As Integer)
               
   ' placer ici la description + signification des paramètres
          
    Dim M As Integer
    
    If N < 2 Then Stop
    If vcible <= T(2) Then
        p1 = 1
    ElseIf vcible >= T(N - 1) Then
        p1 = N - 1
    Else
        p1 = 2
        p2 = N - 1
        Do
            M = (p1 + p2) / 2
            If vcible > T(M) Then p1 = M Else p2 = M
        Loop Until p2 - p1 < 2
    End If
    p2 = p1 + 1
End Sub
Function LocaliseMin(ByVal N As Integer, _
                     ByRef T() As Double, _
            Optional ByRef Tmin As Double) As Integer

               
   ' placer ici la description + signification des paramètres
          
    Dim i As Integer, pos As Integer, v As Double
    
    If N < 1 Then Stop
    Tmin = T(1): pos = 1
    For i = 2 To N
        v = T(i)
        If v < Tmin Then
            pos = i
            Tmin = v
        End If
    Next k
    LocaliseMin = pos
End Function
Function LocaliseMax(ByVal N As Integer, _
                     ByRef T() As Double, _
            Optional ByRef Tmax As Double) As Integer
         
   ' placer ici la description + signification des paramètres
          
    Dim i As Integer, pos As Integer, v As Double
    
    If N < 1 Then Stop
    Tmax = T(1): pos = 1
    For i = 2 To N
        v = T(i)
        If v > Tmax Then
            pos = i
            Tmax = v
        End If
    Next k
    LocaliseMax = pos
End Function


