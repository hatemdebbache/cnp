Attribute VB_Name = "CHERCHE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << CHERCHE >> version 01
'  Requis   :
'  Externes :
'
'  Descrip. : Recherche dans un tableau
'
Option Explicit

Function Vecteur_Max(ByVal k1 As Integer, _
                     ByVal k2 As Integer, _
                     ByRef T() As Double, _
            Optional ByRef p As Integer) As Double

    Dim k As Integer
    Dim m As Double, v As Double
    
    m = -1E+30
    For k = k1 To k2
        v = T(k)
        If v > m Then
            m = v
            p = k
        End If
    Next k
    Vecteur_Max = m
End Function
Function Vecteur_Min(ByVal k1 As Integer, _
                     ByVal k2 As Integer, _
                     ByRef T() As Double, _
            Optional ByRef p As Integer) As Double

    Dim k As Integer
    Dim m As Double, v As Double
    
    m = 1E+30
    For k = k1 To k2
        v = T(k)
        If v < m Then
            m = v
            p = k
        End If
    Next k
    Vecteur_Min = m
End Function
Function Vecteur_MaxAbs(ByVal k1 As Integer, _
                        ByVal k2 As Integer, _
                        ByRef T() As Double, _
               Optional ByRef p As Integer) As Double

    Dim k As Integer
    Dim m As Double, v As Double
    
    m = -1E+30
    For k = k1 To k2
        v = Abs(T(k))
        If v > m Then
            m = v
            p = k
        End If
    Next k
    Vecteur_MaxAbs = m
End Function
Function Vecteur_MinAbs(ByVal k1 As Integer, _
                        ByVal k2 As Integer, _
                        ByRef T() As Double, _
               Optional ByRef p As Integer) As Double

    Dim k As Integer
    Dim m As Double, v As Double
    
    m = 1E+30
    For k = k1 To k2
        v = Abs(T(k))
        If v < m Then
            m = v
            p = k
        End If
    Next k
    Vecteur_MinAbs = m
End Function

